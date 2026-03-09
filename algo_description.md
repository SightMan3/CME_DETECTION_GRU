# Main algorithm pipeline

1. Download 14 hour window of an event
2. Dynamic mask for removing black regions
3. Preprocess images (running diff, to polar, polar map, j-map)
4. Main detection algorithm
5. Process detections (clustering, filtering)
6. Choose best detection
7. Synthetize data based on onset time and velocity approximation.



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

We do things like clipping 1st and 99th percentile, which is just noise. Normalazing the jmap image. Thresholding so that only black and white pixels are present. One of the more important function we use is $\text{mijering}()$. In image processing, the Meijering filter is a specialized algorithm used to detect and enhance continuous, thin, and elongated "ridge-like" or "tubular" structures in an image. Applying it essentially dims the background and highlights objects that look like intersecting lines, webs, or branches. Which is perfect for our usecase and exactly what we need. The last operation that is applied is $\text{closing}()$ in combination with $\text{disk}()$. It creates an artificial disk around a pixel and if there is a pixel in the vicinity of the disk, it connects the two pixels. CMEs on the start and ends of the branch like structure in jamp are sometimes torn and the faintest strucutres are just wandering blobs of pixels, but when connected to the main core structure, it can enhance the approximation of velocity greatly as testing showed.

## 4.2. Detection



