# 2. Data
 
## 2.1 Instrument and observations
 
This study uses white-light coronagraph images from the Large Angle and Spectrometric Coronagraph (LASCO) aboard the Solar and Heliospheric Observatory (SOHO). We use exclusively the C3 detector, whose field of view spans approximately 3.7–32 solar radii. The C3 field of view is the natural choice for kinematic characterisation: it is wide enough that a coronal mass ejection (CME) remains visible over many frames, and far enough from the solar surface that the projected motion of the leading edge is dominated by radial propagation rather than by the initial impulsive acceleration.

Each frame is a 1024 × 1024 image resampled to a common plate scale, with the solar disk centre fixed at pixel (512, 512). Each image sequence is accompanied by instrument metadata retrieved from the Helioviewer archive, from which we take two quantities used throughout the analysis:

- the plate scale $S$, in arcseconds per pixel;
- the apparent solar radius $R_\odot$, in pixels.
All geometric operations are carried out in arcseconds rather than pixels, so that the analysis is insensitive to differences in image scaling between events and epochs.

## 2.2 Event selection and observing window

Events are selected from the SOHO/LASCO CME catalogue maintained by the CDAW Data Center, which provides a catalogued first-appearance time $T_0$ and a linear plane-of-sky speed for each event. For each catalogued event we retrieve a 14-hour image sequence spanning

$$
\left[\,T_0 - 2\,\mathrm{h},\ T_0 + 12\,\mathrm{h}\,\right],
$$

with the lower bound rounded down to the nearest hour. The two-hour lead ensures that the pre-eruption corona is sampled, which is required for the running-difference step (Sect. 3.4) and provides a quiescent baseline against which the eruption is detected. The twelve-hour trailing window is long enough for a CME of even modest speed to traverse the full C3 field of view; a 250 km s$^{-1}$ front crosses the field in roughly six hours.

Each frame carries its acquisition date and time, from which we construct the cumulative time array

$$
\tau_i = t_i - t_0, \qquad i = 0, 1, \ldots, N-1,
$$

used both for temporal regularisation (Sect. 3.5) and for the conversion of detected onset indices back to universal time (Sect. 3.9).
 
## 2.3 Cadence and data gaps
 
The nominal C3 cadence is 12 minutes. Real sequences deviate from this: telemetry dropouts, calibration frames, spacecraft keyhole periods and the reduced cadence of older observations all introduce irregular gaps. This is consequential, because the conversion from a slope measured in image coordinates to a physical velocity requires knowing the temporal width of one pixel along the time axis of the derived time–distance maps. A sequence in which half of the expected frames are missing would yield speeds overestimated by approximately a factor of two if the nominal cadence were assumed uncritically. We therefore quantify the cadence of every sequence explicitly and correct for it where necessary (Sect. 3.5).
 
## 2.4 Validation sample
 
Algorithm performance is assessed on a sample of events spanning a wide range of catalogued speeds (approximately 530–1440 km s$^{-1}$) and angular widths, drawn from several phases of the solar cycle. For each event the CDAW linear speed serves as the reference value. The sample composition, together with the resulting per-event errors, is given in Sect. 4.
 
---
 
# 3. Methods
 
## 3.1 Overview
 
The method estimates two quantities per event: the plane-of-sky propagation speed of the CME leading edge, and the time at which that leading edge entered the C3 field of view. The estimate is obtained not from a single measurement but from an ensemble of independent measurements made at 72 position angles around the corona, which are subsequently associated into candidate events and ranked. The rationale is that a genuine CME produces a spatially coherent signature — consistent onset time and consistent speed across a contiguous range of position angles — whereas noise, cosmic-ray hits and narrow transients such as jets do not.
 
The processing chain comprises six stages:
 
1. removal of instrumental features by geometric masking (Sect. 3.2);
2. transformation of each frame to polar coordinates (Sect. 3.3);
3. construction of running-difference images and time–distance (J-) maps (Sect. 3.4);
4. regularisation of the temporal sampling (Sect. 3.5);
5. enhancement, segmentation and kinematic fitting of propagating fronts at each position angle (Sects. 3.6–3.7);
6. association of individual detections into candidate events, and ranking of those candidates (Sects. 3.8–3.9).
## 3.2 Geometric masking

A raw C3 frame contains four distinct regions, only one of which carries usable coronal signal:

1. the annular field of view containing the corona;
2. the occulting disk, which blocks direct photospheric light;
3. the occulter support pylon, a radial arm crossing the field of view;
4. the black corners of the rectangular detector lying outside the circular field.
Regions 2–4 are not merely uninformative but actively harmful: they are large, contiguous, and have near-zero intensity, and the segmentation stage of the algorithm would otherwise treat their sharp boundaries as legitimate structure. They are therefore removed before any further processing (Sect. 3.2).
 
The pylon poses a particular difficulty because its position angle is not constant. SOHO performs periodic roll manoeuvres, so the pylon appears at different orientations in different events, and its location must be determined per event rather than assumed.
 
For each frame we construct pixel coordinates relative to the disk centre and convert them to arcseconds,
 
$$
x = (i - x_c)\,S, \qquad y = (j - y_c)\,S,
$$
 
and then to polar coordinates,
 
$$
r = \sqrt{x^2 + y^2}, \qquad \theta = \arctan\!\left(\frac{y}{x}\right) \bmod 2\pi ,
$$
 
with $\theta$ wrapped to $[0, 2\pi)$. Working in polar coordinates is what makes the masking tractable: the occulter is a constraint on $r$ alone, the pylon a constraint on $\theta$ alone, and the detector edge again a constraint on $r$ alone. Each of the three masked regions is a coordinate-aligned rectangle in the $(r, \theta)$ plane.
 
Two radial constants are derived once per event and used consistently throughout the analysis:
 
$$
R_{\min} = 4.7\,R_\odot S, \qquad R_{\max} = 512\,S .
$$
 
Masking cculting disk requires all pixels with $r < R_{\min}$ to be masked. The occulter itself subtends 3.7 $R_\odot$; the margin to 4.7 $R_\odot$ accommodates Fresnel diffraction at the occulter edge, which produces a bright diffraction fringe and associated stray light that would otherwise be interpreted as coronal structure.
 
Edge of detector are pixels with $r > R_{\max}$, that are those that lie outside the circular field of view and are masked.
 
Occulter pylon has always fixed width, however its position, specifically the angle, is changing per event. An angular window of width $0.07\pi$ (approximately 12.6°) is swept in steps of $0.01\pi$ through the full circle, restricted radially to the annulus $R_{\min} < r < R_{\max}$ so that the occulter does not contribute. For each window position $\phi$ we compute the fraction of near-zero pixels,
 
$$
D(\phi) = \frac{\left|\{p \in W(\phi) : p < 5\}\right|}{|W(\phi)|},
$$
 
where $W(\phi)$ is the set of pixels in the window. A threshold of 5 rather than 0 is used because image compression leaves the occulted hardware at small but non-zero intensities. The pylon position angle is taken as $\hat{\phi} = \arg\max_\phi D(\phi)$. The window width was fixed empirically and does not vary between events, since the pylon has fixed physical dimensions. The search is performed once per event, on the first frame, and the resulting orientation applied to the whole sequence.
 
The three conditions are combined into a single mask,
 
$$
M = \left\{ r < R_{\min} \right\} \cup \left\{ r > R_{\max} \right\} \cup \left\{ (\hat\phi - 0.07)\pi < \theta < \hat\phi\pi \right\},
$$
 
and masked pixels are set to NaN rather than zero. This distinction matters: zero is a valid intensity and would propagate through subsequent averaging and filtering operations as if it were a measurement, whereas NaN is excluded from all statistics.
 
## 3.3 Polar transformation
 
Each masked frame is resampled onto a regular polar grid of 512 radial and 360 angular samples, spanning $[R_{\min}, R_{\max}]$ in radius and $[0, 2\pi)$ in position angle, giving an angular resolution of exactly 1° per column. For each grid node the corresponding detector coordinates are
 
$$
X = x_c + \frac{r\cos\theta}{S}, \qquad Y = y_c + \frac{r\sin\theta}{S},
$$
 
and the intensity is obtained by bilinear interpolation, with samples falling outside the valid region assigned NaN.
 
The resulting map has radius as its first axis and position angle as its second. The radial grid begins exactly at $R_{\min}$, so the radial pixel scale is
 
$$
\Delta r = \frac{R_{\max} - R_{\min}}{512}
$$
 
by construction. Using the same bounds for the resampling grid and for the subsequent velocity conversion (Sect. 3.7) removes a source of systematic error: if the two are defined inconsistently, every derived speed is biased by their ratio.
 
## 3.4 Running differences and time–distance maps
 
Static structure — the F-corona, stray light, detector artefacts — is removed by differencing consecutive polar maps,
 
$$
\Delta P_i = P_i - P_{i-1}, \qquad i = 1, \ldots, N-1,
$$
 
leaving $N-1$ difference maps. Differences are clipped at zero, retaining only brightenings. The trailing depletion behind a propagating front is thereby discarded; this halves the amount of structure available, but removes the negative counterpart of each front, which otherwise appears in the segmentation stage as a second, spatially offset ridge that competes with the true leading edge.
 
From this sequence we construct time–distance maps, conventionally termed J-maps. For a given position angle $\alpha$, a radial profile is extracted from each difference map and the profiles are stacked chronologically, producing a two-dimensional array
 
$$
J_\alpha(i, k), \qquad i = 0,\ldots,N-2 \ \ (\text{time}), \qquad k = 0,\ldots,511 \ \ (\text{radius}).
$$
 
Rather than extracting a single 1°-wide column, the radial profile is obtained by averaging a 5°-wide angular swath centred on $\alpha$, excluding masked samples from the mean. A single-degree profile is dominated by photon and readout noise, and a swath that intersects a masked region can be lost entirely; averaging over five columns raises the signal-to-noise ratio of the ridge at negligible cost in angular resolution, and matches the 5° sampling interval used for the position-angle scan, so that successive swaths tile the corona without overlap.
 
In a J-map, a radially propagating front appears as an inclined, elongated feature: its position in radius increases monotonically with time. The inclination of this feature encodes the propagation speed, and its extrapolation to the inner edge of the field of view encodes the entry time. A stationary or impulsive brightening, by contrast, appears as a feature of near-zero inclination. Isolating and characterising these inclined features is the core of the detection problem.
 
## 3.5 Temporal regularisation
 
The conversion from J-map inclination to physical velocity presumes a known temporal pixel width. We test this presumption for every sequence by comparing the number of observed sampling intervals against the number expected at nominal cadence over the same span:
 
$$
\kappa = \frac{N - 1}{\left(\tau_{N-1} - \tau_0\right) / \Delta t_{\mathrm{nom}}}, \qquad \Delta t_{\mathrm{nom}} = 720\ \mathrm{s}.
$$
 
A value $\kappa = 1$ indicates uninterrupted nominal sampling; smaller values indicate missing frames.
 
Sequences with $\kappa \geq 0.7$ are analysed as observed. For $\kappa < 0.7$ the J-map is resampled onto a uniform 12-minute grid by linear interpolation along the time axis, using the true frame timestamps as abscissae. Detection then proceeds on the regularised map, and the resulting onset index is rescaled to the original time base by the ratio of map heights before conversion to universal time. Without this rescaling the recovered onset index refers to the interpolated grid and maps to the wrong frame.
 
The threshold $\kappa = 0.7$ is a compromise: interpolation is itself lossy, blurring the very ridges that are subsequently fitted, so it is applied only where the cadence error it corrects exceeds the degradation it introduces.
 
## 3.6 Ridge enhancement and segmentation
 
Each J-map is processed to isolate elongated, ridge-like features. The sequence of operations is as follows.
 
Non-finite values are replaced by zero, and the intensity distribution is clipped at its 1st and 99th percentiles to suppress cosmic-ray hits and defective pixels, then linearly rescaled to $[0, 1]$.
 
Contrast is then equalised using contrast-limited adaptive histogram equalisation (CLAHE), with a clip limit of 0.03. Local rather than global equalisation is essential here, because coronal brightness falls steeply with heliocentric distance: a global stretch is dominated by the bright inner corona and leaves faint outer-corona structure indistinguishable from background. CLAHE applies the stretch within local tiles, so that a faint front at 20 $R_\odot$ receives comparable enhancement to a bright one at 5 $R_\odot$. The clip limit is deliberately conservative; larger values amplify noise into spurious ridge-like texture.
 
Ridge structures are then enhanced using the Meijering neuriteness filter, evaluated at scales $\sigma = 1, 2, 3$ pixels and configured for bright ridges on a dark background. This filter responds to elongated, tubular structures and suppresses both isotropic blobs and extended flat regions, which is well matched to the appearance of a propagating front in a J-map. The scale range corresponds to the observed range of front widths.
 
The filter response is thresholded at its 97th percentile to produce a binary map, and a morphological closing with a disk-shaped structuring element of radius 2 pixels is applied. Closing is required because fronts frequently appear fragmented, particularly at their faint outer extremities, where the ridge breaks into disconnected patches; without closing these are either discarded as undersized or fitted separately, both of which bias the inferred speed. Reconnecting them to the main structure substantially improves speed recovery, at the cost of slightly degrading the recovered spatial extent — an acceptable trade, since speed is the quantity of primary interest for downstream arrival-time modelling.
 
Connected components are then labelled, and components smaller than 40 pixels in area are discarded as noise.
 
## 3.7 Leading-edge fitting and kinematics
 
For each surviving component we determine the inclination of its leading edge. The component is traversed column by column in radius; for each radial column the minimum and maximum time coordinates are recorded, tracing the two temporal boundaries of the feature. A straight line is fitted to each boundary, and the boundary with the smaller absolute inclination $|m|$ — equivalently, the faster of the two — is taken as the leading edge.
 
We use the leading edge rather than the orientation of the component as a whole because the latter is determined by the bulk of the feature, which corresponds to the diffuse body of the ejection, whereas the propagation speed of interest is that of the front.
 
The boundary fits are performed with the RANSAC estimator (residual threshold 1 pixel). Least-squares fitting is unsuitable here because the morphological closing applied in Sect. 3.6 thickens the component irregularly, and a least-squares fit is displaced quadratically by such excursions. RANSAC fits to a maximal consensus subset and excludes boundary pixels lying more than one pixel from the consensus line, recovering the underlying straight edge. Where a component is too small for a consensus set to be established, the fit reverts to ordinary least squares.
 
Two rejection criteria are applied before kinematic conversion. Components spanning fewer than 30 radial pixels are discarded, since a genuine CME traverses a substantial fraction of the C3 field of view. Components with $|m| \leq 0.02$ are discarded as effectively horizontal; such features correspond to brightenings appearing across a wide range of heliocentric distances simultaneously, which is not propagation, and would in any case yield unphysical speeds through the reciprocal in Eq. (1).
 
The inclination is converted to a plane-of-sky speed by
 
$$
v_{\mathrm{POS}} = \left|\frac{1}{m}\right| \cdot \frac{\Delta r}{\Delta t} \cdot C, \tag{1}
$$
 
where $\Delta r$ is the radial pixel scale defined in Sect. 3.3, $\Delta t = 720$ s is the temporal pixel width (guaranteed by Sect. 3.5), and $C = 725$ km arcsec$^{-1}$ converts angular to linear distance at 1 AU. The reciprocal arises because the fitted inclination is a temporal gradient with respect to radius, whereas the speed is the radial gradient with respect to time.
 
The onset index is obtained by extrapolating the fitted leading edge to the inner boundary of the field of view,
 
$$
i_{\mathrm{onset}} = \left\lceil \left| t_c + m\,r_c \right| \right\rceil,
$$
 
where $(r_c, t_c)$ is the component centroid. We emphasise that this quantity is the time at which the leading edge crossed $R_{\min} = 4.7\,R_\odot$, not the eruption time at the solar surface; the two differ by the travel time across the occulted region, which is not observable in C3 and which depends on the acceleration profile below the occulter.
 
Detections with $v_{\mathrm{POS}}$ outside the range 250–2000 km s$^{-1}$ are rejected as physically implausible for a CME.
 
The procedure is repeated at 72 position angles spaced 5° apart, yielding a set of independent detections, each labelled by its position angle, speed and onset index.
 
## 3.8 Association of detections
 
Detections arising from a single eruption should extrapolate to a common onset time regardless of position angle. Detections are therefore sorted by onset index, and a new candidate event is initiated wherever the gap to the preceding detection exceeds one index unit (one frame, nominally 12 minutes). Candidates comprising ten or fewer detections are discarded, corresponding to an angular extent of at most 50°; ejections of interest for terrestrial impact have angular widths of at least 60–70°.
 
Because a single position angle may yield more than one component, each candidate is then reduced so that every position angle contributes exactly one measurement, taking the median onset index and median speed at each angle. The median is preferred to the mean because duplicate detections at a given angle typically consist of one well-fitted track together with a fragment of aberrant inclination.
 
## 3.9 Candidate ranking
 
Each surviving candidate $Q$ is assigned a composite score combining spatial and kinematic coherence,
 
$$
\Theta_Q = 0.40\,\frac{\Delta\theta_Q}{360^\circ} + 0.30\,\frac{F_Q}{100} + 0.15\,\overline{CV}_v + 0.15\,\overline{CV}_t . \tag{2}
$$
 
All four terms are normalised to $[0,1]$, so the weights are directly interpretable. Angular extent carries the largest weight because it is the strongest single discriminant between an ejection and a narrow transient; continuity carries the next largest because a wide but sparsely populated candidate generally reflects chance association rather than a coherent front.
 
**Angular extent.** The angular extent is not the simple difference between extreme position angles, which fails for candidates spanning the 0°/360° discontinuity. Instead, gaps between all adjacent sorted position angles are computed, including the wrap-around gap $(360° - \theta_{\mathrm{last}}) + \theta_{\mathrm{first}}$, and the extent is defined as
 
$$
\Delta\theta_Q = 360° - \max(\text{gaps}),
$$
 
i.e. the full circle less the single largest void. This is exact for both limb events and halo events.
 
**Fill factor.** The fill factor quantifies angular continuity: the fraction of the expected detections within the candidate's angular extent that were actually recovered. With a 5° sampling interval, a candidate of extent $\Delta\theta_Q$ should contain
 
$$
\varepsilon_Q = \frac{\Delta\theta_Q}{5°} + 1
$$
 
detections, the unit term accounting for the inclusive endpoints of the span, and
 
$$
F_Q = \min\left(100,\ \frac{|Q|}{\varepsilon_Q} \times 100\right).
$$
 
**Kinematic coherence.** A physical ejection exhibits comparable speeds and onset times across its angular span. Coherence is measured by the coefficient of variation of each quantity,
 
$$
CV_v = \frac{\sigma_v}{\mu_v}, \qquad CV_t = \frac{\sigma_t}{\mu_t},
$$
 
inverted and bounded so that larger values indicate greater coherence,
 
$$
\overline{CV}_x = \max(0,\ 1 - CV_x).
$$
 
The candidate with the highest $\Theta_Q$ is adopted as the detected event.
 
## 3.10 Final estimates
 
**Speed.** The representative speed of the adopted candidate is obtained by ordering its per-angle speeds by position angle, applying a centred running median of width 3 to suppress isolated outliers, and taking the maximum of the smoothed sequence. The maximum, rather than the mean, is taken because the fastest portion of the front corresponds to the nose of the ejection, which is the component relevant to propagation; slower flank measurements reflect projection of lateral expansion rather than radial motion.
 
**Projection correction.** The speed obtained from Eq. (1) is a plane-of-sky quantity and systematically overestimates the true radial speed for wide events, in which the apparent expansion contains a substantial contribution from lateral rather than radial motion. We apply an empirical correction
 
$$
v = v_{\mathrm{POS}} \left[1 - \eta \min\!\left(\frac{\Delta\theta - 120°}{240°},\, 1\right)\right], \qquad \Delta\theta > 120°, \tag{3}
$$
 
with $\eta = 0.15$; narrower events are left uncorrected, since they lie close to the plane of sky and a correction would introduce rather than remove bias. The correction ramps linearly with angular extent to a maximum reduction of 15% for a full halo. This is deliberately conservative: a rigorous geometric deprojection requires the source region location on the disk, which is not available from coronagraph data alone, and the parametrisation of Eq. (3) should be regarded as a first-order empirical adjustment rather than a physical deprojection.
 
**Onset time.** The onset index of the adopted candidate is taken as the mean over its constituent position angles and mapped to the corresponding frame in the sequence, from which the universal time of field-of-view entry is read directly.
 
## 3.11 Synthetic velocity profiles
 
The recovered speed and onset time serve as initial conditions for a drag-based propagation model, from which a continuous velocity profile between the coronagraph field of view and 1 AU is synthesised. Under the assumption that the dominant force on the ejection is aerodynamic drag against the ambient solar wind, the velocity evolves as
 
$$
v(t) = \frac{v_0 - w}{1 + \gamma (v_0 - w) t} + w,
$$
 
with corresponding radial distance
 
$$
r(t) = \frac{1}{\gamma}\ln\!\left(1 + \gamma(v_0 - w)t\right) + wt + r_0,
$$
 
where $v_0$ is the speed recovered by the detection algorithm, $w$ the ambient solar wind speed, $\gamma$ the drag parameter and $r_0$ the initial heliocentric distance. The transit time is obtained by inverting the second expression at $r = 1$ AU. The drag parameter is modulated by the recovered angular extent, wider (and hence more massive) ejections being less strongly decelerated.
 
Because these profiles are analytic, they are free of the observational noise that affects in-situ measurements, and provide a dense, uniformly sampled feature set for the downstream arrival-time prediction model.
 
## 3.12 Parameter summary
 
| Parameter | Symbol | Value | Section |
|---|---|---|---|
| Observing window | — | $T_0 - 2$ h to $T_0 + 12$ h | 2.2 |
| Inner mask radius | $R_{\min}$ | 4.7 $R_\odot$ | 3.2 |
| Outer mask radius | $R_{\max}$ | 512 px | 3.2 |
| Pylon window width | — | $0.07\pi$ (12.6°) | 3.2 |
| Pylon search step | — | $0.01\pi$ | 3.2 |
| Polar grid | — | 512 radial × 360 angular | 3.3 |
| J-map swath width | — | 5° | 3.4 |
| Position-angle sampling | — | 5° (72 angles) | 3.7 |
| Nominal cadence | $\Delta t_{\mathrm{nom}}$ | 720 s | 3.5 |
| Cadence threshold | $\kappa$ | 0.7 | 3.5 |
| CLAHE clip limit | — | 0.03 | 3.6 |
| Meijering scales | $\sigma$ | 1–3 px | 3.6 |
| Segmentation threshold | — | 97th percentile | 3.6 |
| Closing radius | — | 2 px | 3.6 |
| Minimum component area | — | 40 px | 3.6 |
| Minimum track length | — | 30 px | 3.7 |
| Minimum inclination | $|m|$ | 0.02 | 3.7 |
| Speed bounds | — | 250–2000 km s$^{-1}$ | 3.7 |
| Arcsecond conversion | $C$ | 725 km arcsec$^{-1}$ | 3.7 |
| Minimum detections per candidate | — | > 10 | 3.8 |
| Score weights | — | 0.40 / 0.30 / 0.15 / 0.15 | 3.9 |
| Projection correction cap | $\eta$ | 0.15 | 3.10 |
 
---
 
# 4. Results — outline
 
*(Not derivable from the current material; the following is a proposed structure.)*
 
**4.1 Detection performance.** Fraction of catalogued events for which a candidate was recovered at all, and the number of spurious candidates rejected by the ranking. Worth reporting the distribution of $\Theta$ for accepted versus rejected candidates, since it justifies the weighting in Eq. (2).
 
**4.2 Speed accuracy.** Per-event comparison against the CDAW linear speeds: a table of catalogued versus recovered speed, and a scatter plot with the 1:1 line. Report mean absolute error, but also the signed bias — a systematic over- or under-estimate is diagnostically different from scatter, and the projection correction of Eq. (3) is precisely a bias correction. A residual-versus-angular-width plot will show whether $\eta = 0.15$ is well chosen.
 
**4.3 Onset time accuracy.** Difference between the recovered field-of-view entry time and the catalogued first-appearance time. These are not the same quantity, so a systematic offset is expected rather than an error; the informative number is the scatter about that offset.
 
**4.4 Sensitivity to processing choices.** An ablation is the strongest evidence you can present, and you already have most of the machinery: report MAE with least-squares rather than RANSAC edge fitting, with a single-column rather than 5° swath, with and without the cadence regularisation, and with and without the projection correction. Each of these was adopted on the grounds that it improved results, and the ablation is what converts that claim into evidence.
 
**4.5 Failure modes.** Events where the method fails, and why — overlapping successive eruptions, very slow events that do not clear the field within the window, streamer-blowout events with poorly defined fronts.
 
---
 
## Notes on the transition from technical description to paper
 
Three points that matter for the write-up:
 
1. **The parameters in Sect. 3.12 were tuned on the validation sample of Sect. 2.5.** Reporting an error metric on the same events used to select $\eta$, $\kappa$, the score weights and the thresholds overstates performance. If the sample can be enlarged, split it; if not, state plainly that the quoted errors are in-sample and should be read as a lower bound on the true error.
2. **The onset time is a field-of-view entry time, not an eruption time.** This is stated explicitly in Sect. 3.7 and again in Sect. 3.10, and should not be quietly conflated with the catalogued $T_0$ in the results.
3. **The score weights of Eq. (2) are currently asserted rather than derived.** A sentence on how they were arrived at — grid search, manual tuning against known events, or physical reasoning — is the kind of thing a referee will ask for.