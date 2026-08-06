# Main algorithm pipeline
 
1. Download 14 hour window of an event (+ metadata)
2. Dynamic mask for removing black regions
3. Preprocess images (to polar, running diff, j-map, cadence resampling)
4. Main detection algorithm
5. Process detections (clustering, filtering, scoring, choosing best detection)
6. Synthetize data based on onset time and velocity approximation.
The whole thing from step 2 to step 5 lives inside a single function `run_pipeline(dir_path)`. It takes a directory that contains the PNGs plus an `images_metadata.json` and returns one dictionary describing the best CME candidate found in that directory. Everything is a closure inside that function, because almost every helper needs the per-event metadata (scale, solar radius, pylon angle), and this way we do not have to pass it around everywhere.
 
There are two run modes in the entry point:
- **Normal mode** – optionally wipes `./data_processed/lasco/c3/`, downloads a fresh event, and runs the pipeline once.
- **Test mode** – iterates over every subdirectory of `./data_processed/testing_images_c3/`, runs the pipeline on each one, collects all the results, and finally compares the predicted speeds against the CDAW catalog speeds hardcoded in `real_speeds` and prints the MAE. This is our regression check, it tells us immediately if a change to the algorithm made things better or worse across all test events.
# 1. Download 14 hour window of an event
We start of course by downloading batch of images, for testing and development we take a date in form of YYYY-MM-DD hh:mm:ss in which the event was detected, lets name it $OT$ (Onset Time), its start documented in the SOHO LASCO https://cdaw.gsfc.nasa.gov/CME_list/ catalog. We download and use for detection LASCO c3 images. LASCO is a coronograph on SOHO sun observatory sattelite. Start of the download sequence is calculated as $S = OT - 2h$ (floored to a whole hour) and the end as $E = OT + 12h$. The images themselves are pulled by the `Lasco` extractor class.
 
Along with these images we download also a json file containing metadata regarding the images. This is a separate request to the Helioviewer API endpoint `v2/getClosestImage` with `sourceId=5` (which is LASCO C3), asked for the closest image to the start of our window, and saved as `images_metadata.json` next to the PNGs. From it we use two fields:
- `scale` → arcsec per pixel, how many arcseconds one pixel covers,
- `rsun` → radius of the Sun in pixels.
In this arcsecond scale most of the masking calculations are done.
 
The filenames themselves are also data. They are in the form `YYYYMMDD_HHMM_..._....png`, and we parse them in two places: once to build the cumulative time array of the sequence, and once at the very end to translate a detected onset **index** back into a real date and time.
 
# 2. Dynamic mask for masking out black regions of images.
**LASCO c3 image consits of few main things**:
1. Main circular view.
2. Occulter - disk in front of the sensor to block out the Sun so only solar corona remains.
3. Occulter pylon.
4. Black edges because images are rectangular but view from the coronograph is circular.
The only thing we want to work with is the circular visible region of space with circular corona. We need to replace 0 value pixels with NaN values because the amount of the black pixels would make further algorithms lose to the black pixel because they would use it as valid data.
 
Most important thing to do before we can start masking is creating version of the image in polar cooridnates. We have two matrices $x_{rel}$ and $y_{rel}$. These matrices define x and y coordinates of pixels relative to the center of an image (the center is fixed at pixel 512, 512, because the C3 frames are 1024x1024). We convert them to arcseconds by multiplying by scale.
 
$$
x_{arcsec} = x_{rel}\cdot\text{scale}
$$
$$
y_{arcsec} = y_{rel}\cdot\text{scale}
$$
 
With these matrices we can than use trigonometry to get polar cooridnate matrices:
$$
R_{arcsec} = \sqrt{x_{arcsec}^2 + y_{arcsec}^2} 
$$
$$
\theta_{arcsec} = \tan^{-1}\left(\frac{y_{arcsec}}{x_{arcsec}}\right)
$$
 
So clasicaly we convert cartesian x,y coordinates to polar r, theta coordinates. This will be very helpful to us as we want to accurately in real arcsecond unit define the region of the occulter and the pylon. However currently the $\theta_{arcsec}$ goes from $-\pi$ to $\pi$, which is not really ideal we want to go nicely from $0$ to $2\pi$. That is why we have to do:
 
$$
\theta_{arcsec} \equiv \theta_{arcsec} \mod 2\pi
$$
 
## Two radial constants
Two constants are derived once per event from the metadata and then reused everywhere. They are important because they define the radial extent of everything downstream:
 
$$
R_{min} = 4.7\cdot R_{sun}\cdot\text{scale}, \qquad R_{max} = 512\cdot\text{scale}
$$
 
$R_{min}$ is the outer edge of the occulter mask and $R_{max}$ is the edge of the circular field of view. The same pair is used for the polar unroll grid and for the pixel to arcsecond conversion in the velocity formula, so the geometry stays consistent through the whole pipeline. This is a change from an earlier version, where the unroll used the raw min/max of the radius matrix and the velocity conversion used a different, hardcoded pair (3.7 $R_{sun}$ and half the image width). Those two disagreed with each other and with the mask, which biased the speeds.
 
## Masking the occulter
Occulter is the disk that block Sun. It is always in the center of the image and roughly the same size. We just have to create a circular mask in the center in the image. We will use full 2pi radians so we just need to use radius matrix to set the black pixels, all pixels will be NaNs in the radius 4.7 solar radii.
 
$$
image\left[r_{arcsec} < 4.7\cdot R_{sun}\cdot scale\right]
$$
 
We get the solar radius in pixels, convert it to arcseconds and than multiply by 4.7, so we get radius of 4.7 solar radii. Why is that ? Occulter itself is 3.7 solar radii large, however due to the Fresnel effect, where light bends around the occulter and makes distortions we bump it up to 4.7 juse to be safe.  
 
## Masking the pylon
The pylon is the arm that hold the occulter. It is just a stick however as it goes outwards from the occulter to the camera rear-view, it kind of becomes a cone. In cartesian coordinates it would be really cumbersome to mask this region however this shape in polar cooridantes originates naturally by selecting full length of r and just picking a angle range, so the total opposite of the occulter mask where we took full 360 degrees and just a part of r. The catch is that the pylon can rotate. Due to some correction the SOHO makes the pylon rotates, so to mask it correctly we need to find it first.
 
We will need to find the pylon mask algorithmically. Let 
$$
P = \left(\theta' < X\pi \right) \land \left(\theta' > \left( X - \frac{7}{100}\right)\pi\right)
$$
 
Where X is a fraction of $\pi$. We have to find $X$ so that it matches the pylon in the image. The pylon of the occulter can move and is not the same for every event. Thus the $X$ is decided by following algorithm:
1. Define a donut slice $P$ as a window which would be later the pylon mask, with angle $\phi$ starting at $0.01\pi$.
2. Increment the angle by $\delta = 0.01\pi$, so $\phi = \phi + \delta$, this way we slide the donut window all the way around to $2\pi$.
3. Calculate % number of black pixels in this donut slice, we denote it as $D = \frac{|P < 5|}{|P|}$. Note it is not strictly 0 but under 5, because JPEG/PNG compression leaves the masked hardware at very low but nonzero values.
4. Record the maximum $D$ and store the angle for that $D$. 
5. Than $\phi = X$ and mask for pylon is than: 
$$
\text{image}\left[\left(\theta' < \phi\pi \right) \land \left(\theta' > \left( \phi - \frac{7}{100}\right)\pi\right)\right]
$$
To be clear why it is not a cone slice but a donut slice is that the main view is in donut shape because the occulter is a hole like in a donut and we do not want to count the pixels that are also counted in occulter mask so we ommit it by slicing off the occulter portion of the pylon mask. Concretely the search window is bounded by $R_{min} < r < R_{max}$.
 
Also the fraction $\frac{7}{100}$ was chosen empirically as what seemed to fit the best, the width of the pylon does not change in any way so its okay. $0.07\pi$ is about 12.6 degrees.
 
The search is done **once per event**, on the first image only, and the resulting $\phi$ is then reused for every frame in the sequence.
 
## Black egdes masking
Masking this is also very easy, just like the occulter. We just take all pixels with radius bigger than 512 pixels. Because the radius of the circular view is 512 pixels. So the mask is:
$$
\text{image}\left[r_{arcsec} > 512\cdot \text{scale}\right]
$$
 
## Putting it all together:
$$
\text{image}\left[
r' < 4.7r_{sun} \lor r' > 512 \lor \left[ \left(\theta' < \phi\pi \right) \land \left(\theta' > \left( \phi - \frac{7}{100}\right)\pi\right)\right]
    \right] = \text{NaN}
$$
 
This is applied to every frame by `prepare_image()`, which returns the masked frame together with its radius and angle matrices, so the unroll step does not have to rebuild them.
 
# 3. Image preprocessing
 
After the images are masked properly, the next big step is to unroll them to polar maps, apply running difference and finally make jmaps out of the running difference. We will break it down step by step:
 
## Unrolling to polar map
What we want to achive is to have rectangular image. Yes we currently have a square but the solar corona is wrapped around a circle, the occulter. For later purpouses we need to have an image whose first axis is distance from the Sun and second is the angle. So that the corona flows out of not a circle but a flat line. 
 
This can be achived by some simple trigonometry:
 
$$
X_{\text{pixel}} = X_{\text{center}} + \frac{R_{\text{arcsec}} \cdot \cos(\Theta)}{S_{\text{arcsec/pixel}}}$$
$$Y_{\text{pixel}} = Y_{\text{center}} + \frac{R_{\text{arcsec}} \cdot \sin(\Theta)}{S_{\text{arcsec/pixel}}}$$
 
These are the two core equations that map the polar image (corona wraped around a circle) to where it should be in the XY plane. Basically we are converting a polar coordinate to cartesian coordinate. Lastly we use a package function $\text{map\_coordinates}()$ with linear interpolation (`order=1`) and `cval=NaN`, that uses these $X_{pixel}$ and $Y_{pixel}$ matrices to look for what brightness it should assign from original $\text{image}[]$ data.
 
The sampling grid is built explicitly:
- radius: 512 evenly spaced values from $R_{min}$ to $R_{max}$,
- angle: 360 evenly spaced values from $0$ to $2\pi$ (endpoint excluded), so exactly 1 degree per column.
The result is transposed before returning, so the final polar map has shape **(512 radii, 360 angles)** — rows are distance from the Sun, columns are position angle. The fact that the radial axis starts exactly at $R_{min}$ matters: it means row 0 of the polar map is the first usable pixel outside the occulter, and the pixel-to-arcsecond conversion used later for velocity is exactly $(R_{max} - R_{min})/512$.
 
## Running difference
That is self explanatory. Lets say we have an image $I$ at time $t$, we just do $\Delta I = I_{t} - I_{t-1}$, and we do this for every image, of course not for the first because it has no predecessor. This makes the stationary objects dissapear such as static noise, and moving objects remain and are amplified in a sense. So for $N$ images we end up with $N-1$ difference maps.
 
One detail: the difference is clipped to $[0, 255]$. That means we keep only the brightening (the front of the CME plowing into the previous frame) and throw away the dimming behind it. It halves the amount of structure in the j-map but it also removes the negative ghost trail that used to confuse the ridge filter.
 
## Jmap
This is the pinnacle of the whole image preprocessing. This step is why we needed to unroll the original c3 image to polar maps. Jmap is an image created from slices at a certain angle. Lets say we have 50 images downloaded for our event. In our unrolled polar map we have on one axis the distance from Sun and on the other the angle. Let take angle $a$, for every difference map we take a slice at angle $a$ and stack them chronologically. Doing this we create an image that has on the vertical axis the image index (so if 50 images, 0 to 49) and on horizontal axis distance from sun, that stays the same because we take the slice from polar map with full length.
 
Two things about the current implementation that are not obvious:
 
**It is not a single slice, it is an averaged swath.** Instead of taking column $a$ only, we take a window of `angular_width = 5` columns centered on $a$ and average them with $\text{nanmean}()$ along the angle axis. A one degree slice is extremely noisy, and a single bad column (or a column sitting on the pylon) can destroy the whole track. Averaging 5 degrees together gives a much cleaner ridge for basically free, and 5 degrees is also our angle sampling step so the swaths tile the corona without overlapping. The `nanmean` is wrapped in a warnings suppression because a swath that is entirely inside a mask is all-NaN and numpy complains.
 
**The time axis is reversed.** We iterate the difference maps backwards (`diff_ims[::-1]`), so row 0 of the j-map is the **latest** image and the last row is the earliest. This is why the CME track slopes down-right instead of up-right, and it is why the onset index has to be converted back at the end (see 5.4). It also changes the sign convention in the intercept formula (see 4.2).
 
So the final j-map has shape **(time, radius)**, where the time axis runs backwards.
 
When CME bursts from the sun, without prejudice to generality lets assume the CME is HALO, at every angle of the solar corona we get a moving object, this appears on a jmap as a linear line or a blob that appears linear and elongated with some angle. This way we can capture the moving features in a one static image. The lines or blobs that are created are than used to calculate the onset time as the intercept of a linear function fitted inside the elongated linear object and its slope as the speed. Currently we are speaking only about one angle, but we have to do this through all angles so we can algorithmically decide if what we are detecting is real CME or just a noise, or some other fenomenon like a jet. How we detect these things is discussed in the next chapter.
 
## Cadence check and resampling
This is new and it fixes a whole class of wrong velocities.
 
The velocity conversion assumes that one pixel on the time axis of the j-map is exactly 12 minutes wide, because that is the nominal LASCO C3 cadence. But real events have gaps — SOHO drops out, there are calibration frames, telemetry gaps, and in older data the cadence is different altogether. If half the frames are missing, one pixel is really 24 minutes, and every speed we compute is off by a factor of two.
 
So before detecting on a j-map we check the cadence:
 
$$
\text{cadence score} = \frac{N_{\text{actual intervals}}}{N_{\text{expected intervals}}} = \frac{N_{\text{images}} - 1}{(t_{\text{last}} - t_{\text{first}})/720}
$$
 
The denominator is how many 12 minute steps *should* fit in the observed time span, the numerator is how many we actually have. A score of 1.0 means perfect cadence.
 
- If the score is $\geq 0.7$, the sequence is dense enough and we detect on the j-map as it is.
- If the score is $< 0.7$, the j-map is resampled onto a perfectly even 12 minute grid before detection. We take the real cumulative timestamps (dropping the first one, since the running difference has $N-1$ rows), build a linear interpolator along the time axis with `interp1d(..., fill_value="extrapolate")`, and evaluate it at every 12 minutes from 0 to the end of the window.
After detecting on the resampled j-map, the onset index is in resampled row units, not original row units, so it is scaled back:
 
$$
t_{\text{onset}} \leftarrow t_{\text{onset}} \cdot \frac{H_{\text{original}}}{H_{\text{resampled}}}
$$
 
where $H$ is the number of rows. Without this the onset index points at the wrong frame and the final timestamp is wrong.
 
# 4. Elliptic detection algorithm.
 
The main detection is based on fitting a linear line of the form $y = mx + b$ to some structure in jmap that is created from moving CME. From this linear line we can approximate the speed of the CME and its onset time. 
 
This algorithm has 3 main parts:
1. More preprocessing
2. Detection
3. Calculation of the approximation
## 4.1 More preprocessing.
In this step we alter the image but not in a way we did until now. It will be only filtering of noise and applying functions that form the structures in the image to a state that gives us better chances at detecting even the faintest of CMEs. Things like amplification, morphing pixels, etc...
 
The chain, in order, is:
 
1. **NaN cleanup** – all NaN / +inf / -inf become 0, because the filters downstream cannot handle them.
2. **Percentile clipping** – clip to the 1st and 99th percentile. The extremes are cosmic ray hits and dead pixels, not corona.
3. **Min–max normalization** to $[0,1]$.
4. **CLAHE** – $\text{equalize\_adapthist}()$ with `clip_limit=0.03`. This is contrast limited adaptive histogram equalization: instead of stretching the contrast globally, it does it in local tiles. This matters a lot for us because the corona is dramatically brighter close to the occulter and almost nothing at the field of view edge, so a global stretch would light up the inner corona and leave the outer part black. With CLAHE a faint CME at 20 $R_{sun}$ gets the same contrast treatment as a bright one at 5 $R_{sun}$. The clip limit is deliberately low, at higher values it starts amplifying pure noise into fake ridges.
5. **Meijering ridge filter** – $\text{meijering}(\text{sigmas}=1..3, \text{black\_ridges}=\text{False})$. In image processing, the Meijering filter is a specialized algorithm used to detect and enhance continuous, thin, elongated "ridge-like" or "tubular" structures. Applying it essentially dims the background and highlights objects that look like intersecting lines, webs, or branches. Which is perfect for our usecase, since a CME in a j-map *is* a tubular structure. `black_ridges=False` because our features are bright on dark. The sigma range 1 to 3 sets the widths of ridges we care about, in pixels.
6. **Thresholding** at the 97th percentile of the filtered image, giving a binary map. The percentile is a parameter of the detection function so it can be tuned per experiment, currently 97 everywhere.
7. **Closing with a disk** – $\text{closing}()$ with $\text{disk}(2)$. It creates an artificial disk around a pixel and if there is a pixel in the vicinity of the disk, it connects the two. CMEs at the start and end of the branch-like structure in a jmap are sometimes torn and the faintest structures are just wandering blobs of pixels, but when connected to the main core structure, it can enhance the approximation of velocity greatly as testing showed. It however tampers a bit with the angular width detection and time of onset, however we need to prioritize accurate velocity detection over angular width as this impacts the outcome of our model the most.
8. **Labeling** – $\text{label}()$ turns the binary map into connected components.
## 4.2. Detection
### Mathematical Principles
The core principle of this algorithm is modeling Coronal Mass Ejection (CME) tracks in J-maps by fitting an inclined ellipse to the detected blobs. These blobs are isolated using a Meijering ridge-detection filter. The fitted ellipse is defined in parametric form:$$\frac{(x\cos\theta + y\sin\theta)^2}{(L/2)^2} + \frac{(y\cos\theta - x\sin\theta)^2}{(W/2)^2} = 1$$
Where:
- $\theta$ is the orientation angle of the ellipse relative to the x-axis.
- $L$ is the length of the major axis.
- $W$ is the length of the minor axis.
The ellipse parameters are what `regionprops()` gives us for free, and they are recorded for every detection, but the velocity is **not** derived from the ellipse orientation. The orientation of the fitted ellipse is dominated by the bulk of the blob, and the bulk of the blob is the diffuse body of the CME, not the front. What we actually care about kinematically is the leading edge. So the ellipse is kept as a descriptor and the slope is measured separately.
 
To determine the kinematics of the CME, we calculate its velocity from the leading edge slope:$$v = \left|\frac{1}{m}\right| \cdot \frac{dr}{dt}$$
Where:
- $m$ is the calculated slope of the leading edge.
- $dr$ is the spatial resolution (change in arcseconds per pixel).
- $dt$ is the temporal resolution (change in time per pixel in seconds).
The reciprocal $1/m$ appears because in the j-map the horizontal axis is distance and the vertical axis is time, so the raw slope is $dt/dr$ and we want $dr/dt$.
 
### Region Extraction & Filtering
Algorithm groups the connected binary pixels into distinct objects using $\text{regionprops}()$, these are the blobs we will do calculations with. Then applies a simple filter that discards any blob with a pixel area smaller than 40 to remove small noise artifacts.
 
### Leading edge slope with RANSAC
For every surviving blob we take its pixel coordinates and, for each unique x (radius) column, find the minimum and maximum y (time). That traces the top and bottom edge of the track.
 
Then we fit a line to each edge. This used to be a plain $\text{polyfit}()$ degree 1 fit, and that was a problem: the morphological closing from step 7 makes the blob "fat" and lumpy, and a least squares fit is pulled around by every lump, since a single outlying column influences the fit quadratically. Now both edges are fit with **RANSAC** (`RANSACRegressor`, `residual_threshold=1.0`, `random_state=42`). RANSAC repeatedly fits a line to a random minimal subset and keeps the model with the most inliers, so pixels more than 1 pixel away from the consensus line simply do not vote. In practice it locks onto the true straight edge of the track and ignores the fat left by `closing()`. The random state is fixed so runs are reproducible.
 
If the blob is too small for RANSAC to find a consensus it raises a `ValueError`, and we fall back to plain `polyfit` for that blob.
 
From the two edge slopes we take:
$$
m = \min\left(|m_{\text{min edge}}|, |m_{\text{max edge}}|\right)
$$
The smaller absolute slope is the faster edge (remember, slope here is time over distance), which is the leading front.
 
### Rejection filters
Two cheap filters before we bother with the physics:
- **Track length**: discard if the blob spans less than 30 pixels in x. A real CME crosses a large part of the C3 field of view; a 10 pixel streak is noise.
- **Near-horizontal**: discard if $m \leq 0.02$. A horizontal track in the j-map means something appearing everywhere along the radius at once, which is not propagation, it is a flash, a cosmic ray, or a stray light artifact. It would also produce an absurd velocity through the $1/m$.
### Full velocity calculation
- **Spatial resolution** ($dr$): the field of view bounds distributed over the 512 radial pixels of the polar map. These are the *same* $R_{min}$, $R_{max}$ used to build the polar grid, which is the fix mentioned in section 2:
$$dr = \frac{R_{max} - R_{min}}{512}$$
- **Temporal resolution** ($dt$): the standard cadence of 12 minutes, converted to seconds:
$$dt = 12 \cdot 60 = 720 \text{ s}$$
This is valid because of the cadence check in section 3 — either the sequence really is at 12 minutes, or it was resampled to be.
- **Final Calculation**: velocity in arcseconds per second, then converted to km/s with `ARCSEC_TO_KM = 725`:
$$v_{\text{arcsec/s}} = \left| \frac{1}{m} \right| \cdot \frac{dr}{dt}$$
$$v_{\text{km/s}} = v_{\text{arcsec/s}} \cdot \text{ARCSEC\_TO\_KM}$$
To filter out unrealistic detections, only blobs with a final velocity between **250 km/s and 2000 km/s** are retained.
 
### Onset index
The estimated onset time index is the intercept of the leading edge line, extrapolated back to $r = 0$:
$$t_{\text{onset}} = \left\lceil |y_0 + m \cdot x_0| \right\rceil$$
 
The **plus** sign is not a typo. Because the j-map time axis is reversed (section 3), a CME moving outward has a genuinely negative slope, but we stored $m$ as an absolute value. Substituting $-m$ into the usual intercept $y_0 - m x_0$ gives $y_0 + m x_0$. So this index is counted from the *end* of the sequence backwards, which is undone in section 5.4.
 
If the blob passes everything, its centroid, axes, orientation, slope, velocity and onset index are appended to the detection list.
 
# 5. Detection processing
The process of detection is applied to only one j-map which we know is a representation of only one angle from the corona. We have to do this for every angle. However applying detection 360x is time consuming. For simplicity we use a 5 degree sample rate, this means we are making only one j-map per 5 degrees, which is 360/5 = 72 j-maps and also roughly 72 detections. Roughly, because we can get no detection for certain angles, or multiple detections for a single angle. And that is where detection processing comes in. We need to know which detections are the CME we are trying to detect and which are noise or some other type of object like a jet.
 
All detections from all angles are collected into one dataframe with an added `angle` column.
 
## 5.1 Clustering
We create clusters of detections based on their approximated onset index. The detections are sorted by `t_onset_idx`, the consecutive difference is taken, and a new cluster starts whenever the gap to the previous detection is larger than 1 index step:
 
$$
\text{cluster id} = \text{cumsum}\left(\Delta t_{\text{onset}} > 1\right)
$$
 
One index step is one j-map row, so about 12 minutes. Detections belonging to the same eruption should all extrapolate back to within a frame or two of each other, regardless of which angle they came from.
 
Then clusters with 10 or fewer detections are dropped, because 10 angles or less cannot be a valid CME. CMEs we want to detect have angular width of at least 60-70 degrees, and even then there is not a big likelihood of a CME with 60 degree AW heading towards Earth.
 
## 5.2 Flattening
A single angle can produce several blobs, and several blobs can land in the same cluster. Before scoring, each cluster is collapsed so that **one angle contributes exactly one row**, by grouping on `angle` and taking the **median** of `t_onset_idx` and of `velocity_km_s`. Median rather than mean, because when an angle produces a duplicate detection it is usually one good track plus one fragment with a wild velocity, and the median ignores it. This also makes the fill factor meaningful, since it needs unique angles to count.
 
## 5.3 Quality score & fill factor
After individual detections are grouped and flattened, the algorithm must determine which cluster represents the most probable CME. This is achieved by calculating a composite Quality Score ($\Theta$) for each cluster. The evaluation relies on assessing the physical structure (Angular Width and Fill Factor) and the kinematic coherence (variance in velocity and onset time) of the cluster.
 
Let $Q$ be a cluster, $Q = \{q_1, q_2, \ldots, q_n\}$ where $q_i$ is a detection at a unique angle. The score is:
 
$$
\Theta_{Q} = 0.4\cdot\frac{\Delta\theta_Q}{360} + 0.3\cdot\frac{FF_Q}{100} + 0.15\cdot\overline{CV}_{Q_v} + 0.15\cdot\overline{CV}_{Q_t}
$$
 
All four terms are normalized to $[0,1]$ so the weights actually mean what they look like. Angular width dominates (0.4) because it is the strongest single discriminator between a CME and a jet, fill factor is next (0.3) because a wide but sparse cluster is usually an accident of clustering, and the two coherence terms split the remaining 0.3 evenly.
 
### True angular width $\Delta \theta$
Instead of simply taking the difference between the maximum and minimum angles (which fails when a CME crosses the 360°/0° boundary), the algorithm calculates the gaps between all adjacent, sorted angles. It also calculates the "wrap-around" gap to account for the circular geometry of the coronagraph:
$$\text{Wrap Gap} = (360^\circ - \theta_{\text{last}}) + \theta_{\text{first}}$$
The true Angular Width is defined by subtracting the single largest empty gap from the full circle:
$$\Delta \theta = 360^\circ - \text{max}(\text{gaps})$$
 
### Fill Factor
**Fill Factor** tells us how continuous the cluster is. When the angular width $\Delta\theta_Q$ of a cluster is for example 340 degrees, then we expect roughly $340/5$ detections in it, because 5 degrees is our angular resolution. The expected count is
$$\epsilon_Q = \frac{\Delta\theta_Q}{5} + 1$$
The $+1$ is a fencepost correction: a span covering angles 0, 5, 10 has a width of 10 degrees but contains 3 samples, not 2. Without it, a perfectly filled cluster scores above 100%.
 
$$
FF_Q = \min\left(100,\ \frac{|Q|}{\epsilon_Q} \cdot 100\right)
$$
 
capped at 100 so a slight overcount cannot inflate the score.
 
### Kinematic Coherence: Coefficient of Variation (CV)
A physical CME should exhibit relatively consistent velocities and onset times across its angular span. To measure this coherence, the algorithm calculates the Coefficient of Variation for both velocity ($v$) and onset time ($t$). CV is the ratio of the standard deviation to the mean:
$$CV_v = \frac{\sigma_v}{\mu_v}, \quad CV_t = \frac{\sigma_t}{\mu_t}$$
A lower CV indicates higher coherence (less noise). To translate this into a scoring metric where "higher is better", the CV is inverted and bounded at zero:
$$\overline{CV}_v = \text{max}(0, 1 - CV_v)$$
$$\overline{CV}_t = \text{max}(0, 1 - CV_t)$$
If the mean is zero the CV defaults to 1.0, so the term contributes nothing rather than dividing by zero.
 
These are computed on the raw per-angle velocities, not the smoothed ones, because the point of the term is precisely to measure the scatter.
 
### Representative Velocity Estimation
To report a single, realistic velocity for the entire CME cluster, the algorithm sorts the detections by angle and applies a rolling median filter (window size 3, centered, `min_periods=1`). This smooths out local outliers and noise. The maximum value of these smoothed velocities is taken as the plane-of-sky CME velocity, effectively identifying the fastest moving part of the leading edge (the "nose" of the CME).
 
### Projection effect correction
The velocity we measure is a plane-of-sky velocity. For a CME coming straight at us (which is exactly the case we care about for arrival prediction), the plane-of-sky speed is a projection and it systematically overestimates the true radial speed, because what we see expanding sideways is partly the lateral expansion of a halo, not radial propagation.
 
So the plane-of-sky speed is reduced for wide events:
 
$$
v_{\text{corrected}} = v_{\text{POS}} \cdot \left(1 - 0.15\cdot\min\left(\frac{\Delta\theta - 120}{240}, 1\right)\right), \qquad \Delta\theta > 120°
$$
 
Below 120 degrees the velocity is left untouched — a narrow CME is nearly in the plane of sky already, and correcting it would be wrong. Above 120 degrees the penalty ramps linearly up to a maximum of **15% reduction** for a full 360 degree halo. The 15% cap is empirical, tuned against the catalog speeds in the test harness; it is deliberately conservative because a real geometric deprojection needs the source location on the disk, which we do not have.
 
The reported `velocity_km_s` is this corrected value. The score itself uses the uncorrected spread.
 
## 5.4 From index to timestamp
The onset index is an index into the reversed j-map, so it has to be flipped back and mapped to a real file:
 
$$
ot = \lfloor \mu_t \rfloor, \qquad iot = N_{\text{images}} - ot
$$
 
$iot$ is then used to index into the sorted list of PNG filenames, and the filename is split on `_` to recover the date (`YYYYMMDD`) and time (`HHMM`) of the frame where the CME started. Both indices and both string fields are reported, mainly so the mistake is visible if the flip is ever off by one.
 
## 5.5 Output and best cluster selection
Each cluster becomes a dictionary entry with: the two CVs, the corrected velocity, the onset index and its inverse, the onset date and datetime strings, the angular width, the fill factor, and $\Theta$.
 
If no clusters survived at all, a single placeholder entry is returned with all fields zeroed and $\Theta = -1$, so the caller gets a well-formed dictionary instead of an exception and the test harness can keep going through the remaining events.
 
Naturally the cluster with the highest quality score is considered the best candidate for the CME we are detecting, and its approximated velocity and onset time will be used for further calculations and considered as real. `run_pipeline` returns only that one dictionary.
 
## 5.6 Notes on the current implementation
A few things worth being aware of when reading the code:
 
- `estimate_angular_width()` computes an angular width directly from the difference image cube (nanmax intensity profile over all frames, 95th percentile threshold, largest-gap logic, with the rule that a gap wider than 90 degrees means a partial CME). It is called in the scoring loop but its result is currently **not used** — both the reported angular width and the projection correction use the fill-factor-derived $\Delta\theta$ instead. It is kept as an alternative estimator to compare against.
- `image_masking()` is superseded by `prepare_image()` and is no longer called.
- The two functions build the angle matrix with different conventions (`prepare_image` does not apply the $\bmod 2\pi$ and flips the sign of the relative y axis compared to the pylon search). Since the pylon angle $\phi$ is found in one convention and applied in the other, the pylon mask may be landing at the mirrored position rather than on the pylon. This is worth verifying by overlaying the mask on a frame.
- The cadence score is recomputed inside the angle loop although it does not depend on the angle, so it is evaluated 72 times per event instead of once.
# 6. Synthetic velocity
After approximating the velocity and onset time of the CME we can use these values to calculate the expected evolution of the CME's velocity over time until it gets to Earth. This part lives outside `elliptic_cactus.py`, in the data preparation module.
 
The main idea is to create a time series of angular width and of velocity. The velocity is modeled with a drag based model. We use the initial velocity calculated with the CME detection algorithm as initial condition:
 
$$
v(t) = \frac{v_0 - w}{1 +\gamma(v_0 - w)t} + w
$$
 
where $v_0$ is the initial velocity, $w$ is the wind velocity, $\gamma$ is the drag coefficient which we set to $0.5 \times 10^{-7}$, later adjusted based on angular width: if the CME is halo, that means it has more mass, so it will go through space easier and should experience lower drag compared to CMEs with angular width less than 120 degrees. $t$ is the time elapsed since launch in seconds.
 
The velocity function gives us velocity at time $t$, so we need a time span in which the CME is propagating through space towards Earth. That we calculate by solving this equation for $t$:
$$
r(t) = \frac{1}{\gamma} \ln{\left(1 + \gamma(v_0 - w)t\right)} + tw + r_0.
$$
 
These attributes will then be used for training prediction models. Because this velocity attribute is purely synthetic there is no noise, which maximizes the information yield for the prediction model.