# Main algorithm pipeline

1. Download 14 hour window of an event
2. Dynamic mask for removing black regions
3. Preprocess images (running diff, to polar, polar map, j-map)
4. Main detection algorithm
5. Process detections (clustering, filtering, choosing best detection)
6. Synthetize data based on onset time and velocity approximation.



# 1. Download 14 window of an event
We start of course by downloading batch of images, for testing and development we take a date in form of YYYY-MM-DD hh:mm:ss in which the event was detected, lets name it $OT$ (Onset Time), its start documented in the SOHO LASCO https://cdaw.gsfc.nasa.gov/CME_list/ catalog. Data is downloaded with hapi client. We download and use for detection LASCO c3 images. LASCO is a coronograph on SOHO sun observatory sattelite. Start of the download sequence is calculated as $S = OT - 2h$ and the end as $E = OT + 12h$. 

Along with these images is downloaded also a json file containing metadata regarding the imges. Information like distance from the Sun, diamater of the Sun. Image height, width and scale which helps us convert pixels to arcsec per pixel which says how many degrees in arcseconds is one pixel. In this scale most of the masking calculations are done.

# 2. Dynamic mask for masking out black regions of images.
**LASCO c3 image consits of few main things**:
1. Main circular view.
2. Occulter - disk in front of the sensor to block out the Sun so only solar corona remains.
3. Occulter pylon.
4. Black edges because images are rectangular but view from the coronograph is circular.

The only thing we want to work with is the circular visible region of space with circular corona. We need to replace 0 value pixels with NaN values because the amount of the black pixels would make further algorithms lose to the black pixel because they would use it as valid data.

Most important thing to do before we can start masking is creating version of the image in polar cooridnates. We have two matrices $x_{rel}$ and $y_{rel}$. These matrices define x and y coordinates of pixels relative to the center of an image. We convert them to arcseconds by multiplying by scale.

$$
x_{arcsec} = x_{rel}\cdot\text{scale}
$$
$$
y_{arcsec} = y_{rel}\cdot\text{scale}
$$

With these mmatrices we can than use trigonometry to get polar cooridnate matrices:
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

## Masking the occulter
Occulter is the disk that block Sun. It is always in the center of the image and roughly the same size. We just have to create a circular mask in the center in the image. We will use full 2pi radians so we just need to use radius matrix to set the black pixels, all pixels will be NaNs in the radius 4.7 solar radii.

$$
image\left[r_{arcsec} < 4.7\cdot R_{sun}\cdot scale\right]
$$

We get the solar radius in pixels, convert it to arcseconds and than multiply by 4.7, so we get radius of 4.7 solar radii. Why is that ? Occulter itself is 3.7 solar radii large, however due to the Fresnel effect, where light bends around the occulter and makes distortions we bump it up to 4.7 juse to be safe.  

## Masking the pylon
The pylon is the arm that hold the occulter. It is just a stick however as it goes outwards from the occulter to the camera rear-view, it kind of becomes a cone. In cartesian coordinates it would be really cumbersome to maks this regions however this shape in polar cooridantes originates naturally by selecting full length of r and just picking a angle range, so the total opposite of the occulter mask wehre we took full 360 degrees and just a part of r. The catch is that the pylon can rotate. Due to some correction the SOHO makes the pylon rotates, so to mask it correctly we need to find it first.
e will need to find the pylonn mask algorithmically. Let 
$$
P = \left(\theta' > X\pi \right) \land \left(\theta' < \left( X - \frac{7}{100}\right)\pi\right)
$$

Where X is a fraction. Angle in radians is $\frac{a}{b}\pi$. We have to find the fraction $X = \frac{a}{b}$ so that it matches the pylon in the image. The pylon of the occulter can move and is not the same for every event. Thus the $X$ is decided by following algorithm:
1. Define a donut slice $P$ as a window which would be later the pylon mask with angle $\phi$ starting at 0 radians.
2. Increment the angle by some small $\delta$, so $\phi = \phi + \delta$, this way we slide the donut window
3. caculate % number of black pixels in this donut slice, we denote it as $D = \frac{|P < 5|}{|P|}$.
4. Record the maximum $D$ and store the angle for that $D$. 
5. Than $\phi = X$ and mask for pylon is than: 
$$
\text{image}\left[\left(\theta' > \phi\pi \right) \land \left(\theta' < \left( \phi - \frac{7}{100}\right)\pi\right)\right]
$$

To be clear why it is not a cone slice but a donut slice is that the main view is in donut shape because the occulter is a hole like in a donut and we do not want to count the pixels that are also counted in occulter mask so we ommit it by slicing off the occulter portion of the pylon mask.

Also the fraction $\frac{7}{100}$ was chosen empirically as what seemed to fit the best, the width of the pylon does not change in any way so its okay.

## Black egdes masking
Masking this is also very easy, just like the occulter. We just take all pixels with radius bigger than 512. Because the radius of the circular view is 512 pixels. So the mask is:
$$
\text{image}\left[r_{arcsec} > 512\cdot \text{scale}\right]
$$

## Putting it all together:
$$
\text{image}\left[

r' < 4.7r_{sun} \lor r' > 512 \lor \left[ \left(\theta' > \phi \right) \land \left(\theta' < \left( \phi - \frac{7}{100}\right)\pi\right)\right]

    \right]
$$

# 3. Image preprocessing

After the images are masked properly, the next big step is to unroll them to polar maps, apply running difference and finally make jmaps out of the running difference. We will break it down step by step:

## Unrolling to polar map
What we want to achive is to have rectangular image. Yes we currently have a square but the solar corona is wrapped around a circle, the occulter. For later purpouses we need to have an image whch first axis is radius and second id distance from sun. So that the corona flows out of not a circle but a flat line. 

This can be achived som simple trigonometry:

$$
X_{\text{pixel}} = X_{\text{center}} + \frac{R_{\text{arcsec}} \cdot \cos(\Theta)}{S_{\text{arcsec/pixel}}}$$
$$Y_{\text{pixel}} = Y_{\text{center}} + \frac{R_{\text{arcsec}} \cdot \sin(\Theta)}{S_{\text{arcsec/pixel}}}$$

These are the two core equations that map the polar image (corona wraped around a circle) to whre it should be in the XY plane. Basically we are converting a polar coordinate to cartesian coordinate, the polar coordinates are the matrices we used before to create masks. Lastly we use a package function $\text{map\_coordinates}()$ that uses these $X_{pixel}$ and $Y_{pixel}$ matrices to look for what brightness it should assign from original $\text{image}[]$ data.

## Running difference
That is self explanatory. Lets say we have and image $I$ at time $t$, we just do $\Delta I = I_{t-1} - I_t$, and we do this for every image, of course not for the first because it has no predecessor. This makes the stationary object dissapear such as static noise, and moving object remain and are amplified in a sense.

## Jmap
This is the pinnacle of the whole image preprocessing. This step is why we needed to unroll the original c3 image to polar maps. Jmap is an image created from slices at a crtain angle. Lets say we have 50 images downloaded for our event. In our unrolled polar map we have on x-axis degrees from 0 to 360 and on y-axis we have the distance from sun. Let take angle $a$, for every image we take a slice at angle $a$ and put it chronologically in a jmap. So for the first image we take a slice, than second, thrid... etc and stack side by side. Doing this we create an image that has on y-axis number of images, so if 50 than we y-axsi goes from 0-50 and on x-axis we have distance from sun, that stays the same because we take the slice from polar map with full length.

When CME bursts from the sun,  without prejudice to generality lets assume the CME is HALO, at every angle of the solar corona we get a moving object, this appears on a jmap as a linear line or a blob that appears linear and elongated with some angle headed upwards. This way we can capture the moving features in a one static image. The lines or blobs that are created are than used to calculate the onset time as y-intercept of linear function fitted inside the elongated linear object and its slope as the speed. Currently we are speaking only about one angle, but we have to do this through all angles so we can algorithmically decide if what we are detecting is real CME or just a noise, or some other fenomenon like a jet. How we detect these things is discussed in the next chapter.

# 4. Elliptic detection algorithm.

The main detection is based on fitting a linear line of the form $y = mx + b$ to some structure in jamp that is created from moving CME. from this linear line we can approximate the speed of the CME and its onset time. 

This algorithm has 3 main parts:
1. More preprocessing
2. Detection
3. Calculation of the approximation

And also few things about saving the detection in some form but that is not that important.

## 4.1 More preprocessing.
In this step we alter the image but not in a way we did until now. It will be only filtering of noise and applying functions that form the structures in the image to a state that gives us better chances at detecting even the faintest of CMEs. Things like amplification, morphing pixels, etc...

We do things like clipping 1st and 99th percentile, which is just noise. Normalazing the jmap image. Thresholding so that only black and white pixels are present. One of the more important function we use is $\text{mijering}()$. In image processing, the Meijering filter is a specialized algorithm used to detect and enhance continuous, thin, and elongated "ridge-like" or "tubular" structures in an image. Applying it essentially dims the background and highlights objects that look like intersecting lines, webs, or branches. Which is perfect for our usecase and exactly what we need. The last operation that is applied is $\text{closing}()$ in combination with $\text{disk}()$. It creates an artificial disk around a pixel and if there is a pixel in the vicinity of the disk, it connects the two pixels. CMEs on the start and ends of the branch like structure in jamp are sometimes torn and the faintest strucutres are just wandering blobs of pixels, but when connected to the main core structure, it can enhance the approximation of velocity greatly as testing showed. It however tampered little bit with the angular width detection and time of onset in C3, however we need to prioritze accurate velocity detection over angular width as this impacts the outcome of our model the most.

## 4.2. Detection
### Mathematical Principles
The core principle of this algorithm is modeling Coronal Mass Ejection (CME) tracks in J-maps by fitting an inclined ellipse to the detected blobs. These blobs are isolated using a Meijering ridge-detection filter.The fitted ellipse is defined in parametric form:$$\frac{(x\cos\theta + y\sin\theta)^2}{(L/2)^2} + \frac{(y\cos\theta - x\sin\theta)^2}{(W/2)^2} = 1$$
Where:
- $\theta$ is the orientation angle of the ellipse relative to the x-axis.
- $L$ is the length of the major axis.
- $W$ is the length of the minor axis.

To determine the kinematics of the CME, we calculate its velocity. Instead of relying simply on the angle of the major axis ($\theta$), the algorithm calculates the slope of the leading edge of the blob to better represent the propagating front of the CME. The kinematic velocity is then calculated as:$$v = \left|\frac{1}{m}\right| \cdot \frac{dr}{dt}$$
Where:
- $m$ is the calculated slope of the leading edge.
- $dr$ is the spatial resolution (change in arcseconds per pixel).
- $dt$ is the temporal resolution (change in time per pixel in seconds).

To filter out noise and unrealistic detections (such as nearly vertical or horizontal slopes), the algorithm checks if the calculated velocity falls within realistic bounds for a CME. Only blobs with a final velocity between 250 km/s and 2000 km/s are retained.Additionally, the estimated onset time index is derived linearly from the centroid coordinates ($x_0, y_0$) and the leading edge slope:$$t_{\text{onset}} = |y_0 + m \cdot x_0|$$

### Region Extraction & Filtering
Algorithm groups the connected binary pixels into distinct objects using $\text{regionprops}()$ these are the blobs we will do calculations with. Than applies simple filter that discards any blobs with a pixel area smaller than 40 to remove small noise artifacts.

### Main pipeline
For every valid blob, the algorithm extracts its centroid, major/minor axes, and orientation. It then performs the following calculations. Leading Edge Slope Calculation: Iterates through the unique x-coordinates of the blob. It finds the minimum and maximum y-coordinates for each x-slice, effectively tracing the top and bottom edges of the track. It fits a linear polynomial (degree 1) to both edges and selects the slope with the smallest absolute value to represent the leading edge ($m$). Discards tracks shorter than 30 pixels in the x-direction. Discards near-horizontal tracks (slope $\le$ 0.02). Converts the pixel-based slope into real-world units (arcseconds per second) using predefined constants ($R_{arcsec}$ bounds, and an assumed $dt$ of 12 minutes), and subsequently converts it to km/s. 12 minutes because 1 pixel, on the time axis, is 12 minutes wide, and we need to know that for the conversion calculation.
Here is the full calculation of speed in code:
- **Spatial resolution** ($dr$): Calculated based on the field of view bounds ($R_{\text{max}}$ and $R_{\text{min}}$ in arcseconds) distributed over the 512-pixel height of the J-map slice:
$$dr = \frac{R_{\text{max}} - R_{\text{min}}}{512}$$
- **Temporal resolution** ($dt$): Assumes a standard image cadence of 12 minutes, converted to seconds:
$$dt = 12 \cdot 60 = 720 \text{ s}$$
- **Final Calculation**: The velocity is first calculated in arcseconds per second using the leading edge slope ($m$) with the formula we introduced earlier, and then converted to kilometers per second using a predefined conversion constant (ARCSEC_TO_KM):$$v_{\text{arcsec/s}} = \left| \frac{1}{m} \right| \cdot \frac{dr}{dt}$$
$$v_{\text{km/s}} = v_{\text{arcsec/s}} \cdot \text{ARCSEC\_TO\_KM}$$

If the blob's calculated velocity falls within the defined physical bounds (250 km/s to 2000 km/s), the properties of the CME (centroid, axes, slope, velocity, and onset index) are saved and returned in a dictionary.

# 5. Detection processing
The process of detection is applied only one j-map which we know is representation of only one angle from the corona. We have to do this for every angle j-map. However applying detection 360x is time consuming if we have to do over and over. For simplicity we used 5 angle sample rate, this means we are making only one j-map per 5 angles, which is 360/5 = 72 j-maps and also roughly 72 detections. Roughly because we cant get no detections for certain angles, or multiple detections for certain angle. And that's where detection processing comes in. We need to know which detections are the CME we are trying to detect and which are noise or some other types of objects like jets.

## 5.1 Clustering
We create clusters of detections based on thir approximated $OT$. For every detection in a list we calculate to which cluster it falls into. Its just a number 0,1,2,3,4,... It filters out clusters with less than 10 detections, because 10 angles or less cant be a valid cme. CMEs we want to detect have angular width of at least 60-70 degrees, even then there is not big likelyhood of a CME with 60 degree PA heading towards earth.

## 5.2 Quality score & fill factor
After individual detections are grouped into clusters based on their onset times, the algorithm must determine which cluster represents the most probable CME. This is achieved by calculating a composite Quality Score ($\Theta$) for each cluster.The evaluation relies on assessing the physical structure (Angular Width and Fill Factor) and kinematic coherence (variance in velocity and onset time) of the cluster.

Let $Q$ be a cluster, $Q = \{q_1, q_2, \ldots, q_n\}$, where $q_i$ is a detection. We than calculate the score $\Theta$ for the cluster $Q$ as:
$$
\Theta_{Q} = w_1\Delta\theta_Q + w_2CV_{Q_t} + w_3CV_{Q_v}  
$$

**Fill Factor** tells us how continuous is the cluster. When our true angular width $\Delta\theta_Q$ of the cluster $Q$ is for example 340 degrees, then we expect at least $340/5$ (divided by 5 because that is our resolution) detections in the cluster $Q$, we will call this value as expected number of detections $\epsilon_Q = \Delta\theta_Q/5$ and the true number of **unique** detections will be just $|Q|$ the size of the set. 
$$
FF_Q = \frac{|Q|}{\epsilon_Q} \cdot 100
$$
We can incorporate this into our score:
$$
\Theta_{Q} = w_1\Delta\theta_Q + w_4FF_Q + w_2CV_{Q_t} + w_3CV_{Q_v}
$$

Where CV is coeffcient of variation of the clusters velocity and onset time. We invert its value to make a higher value better. Like so:
$$
CV_{Q_t} = \text{max}(0, 1 - CV_{Q_t} )
$$
$$
CV_{Q_v} = \text{max}(0, 1 - CV_{Q_v} )
$$

### True angular width $\Delta \theta$
Instead of simply taking the difference between the maximum and minimum angles (which fails when a CME crosses the 360°/0° boundary), the algorithm calculates the gaps between all adjacent, sorted angles. It also calculates the "wrap-around" gap to account for the circular geometry of the coronagraph:
$$\text{Wrap Gap} = (360^\circ - \theta_{\text{last}}) + \theta_{\text{first}}$$
The true Angular Width is defined by subtracting the single largest empty gap from the full 360-degree circle:
$$\Delta \theta = 360^\circ - \text{max}(\text{gaps})$$

### Kinematic Coherence: Coefficient of Variation (CV)
A physical CME should exhibit relatively consistent velocities and onset times across its angular span. To measure this coherence, the algorithm calculates the Coefficient of Variation (CV) for both velocity ($v$) and onset time ($t$). CV is the ratio of the standard deviation to the mean:$$CV_v = \frac{\sigma_v}{\mu_v}, \quad CV_t = \frac{\sigma_t}{\mu_t}$$A lower CV indicates higher coherence (less noise). To translate this into a scoring metric where "higher is better," the CV is inverted and bounded at zero:
$$\text{Norm\_CV}_v = \text{max}(0, 1 - CV_v)$$
$$\text{Norm\_CV}_t = \text{max}(0, 1 - CV_t)$$

### Representative Velocity Estimation

To report a single, realistic velocity for the entire CME cluster, the algorithm sorts the detections by angle and applies a rolling median filter (with a window size of 3). This smooths out local outliers and noise. The maximum value of these smoothed velocities is taken as the true CME velocity, effectively identifying the fastest moving part of the leading edge (the "nose" of the CME).


### Conclusion
Naturally the cluster with highest quality score is considered as the best candidate for the CME we are detecting, and its approximated velocity and onset time will be used for further calculations and will be considered as real.

# 6. Synthetic velocity
After appriximating the velocity and onset time of CME we can use these values to calculate expected evolution of CMEs velocity over time until it gets to earth.

The main idea is to create time series of angular width and of velocity. The velocity will be modeled with drag based model. We will use the intial velocity calculated with the cme detection algorithm as initial condition for the equation for velocity over time:

$$
v(t) = \frac{v_0 - w}{1 +\gamma(v_0 - w)t} + w
$$

where $v_0$ is the initial velocity, $w$ is the wind velocity, $\gamma$ is the drag coefficient which we will set to $0.5 \times 10^{-7}$, this will we be later adjusted based on the angular width of CMEs, if the CME is halo, that means it has more mass, so it will go through space easier and should experience lower drag compared to CMEs with angular width less than 120 degrees. $t$ is the time elapsed since launch in seconds.

The velocity function give us velocity at time t, so we will need a time span in which the cme is propagating through space towards earth. That we will calculate by solving this equation for $t$:
$$
r(t) = \frac{1}{\gamma} \ln{\left(1 + \gamma(v_0 - w)t\right)} + tw + r_0.
$$


These attributes will than further be used for training prediction models. Because this velocity attribute is purely synthetic there is no noise which maximazes the information yield for the prediction model.










