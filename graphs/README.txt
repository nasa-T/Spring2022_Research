Documentation for Graphs:

All graphs were created with the current version of each Python file 
(no further changes to parameters in the Python files have been made following the creation of these graphs — besides commenting out 'plt.savefig(...)')

Note: some graphs will have to change formatting for readability
'~' signifies plots that will change on each run due to randomness inherent to generating the plots


File name; python file of origin; Description

ALMA_WISE_correlation.png; ALMA_WISE_correlation.py; bar graphs displaying potential relationships with ALMA detections and WISE band detections

J-H_H-K_&_WISE_diskClassification.png; disk_fraction.py; representation of distribution of objects on color-color diagrams (J-H_H-K, W1-W2_W2-W3, W1-W2_W3-W4) along with the relevant lines for disk classifications. The JHK diagram uses the main sequence curve from Pecaut & Mamajek (2013), the Interstellar reddening vector from Cardelli et al. (1989), and the Classical T Tauri Star locus from Meyer et al. (1997). W123 diagram possesses lines signifying cutoffs for objects with and without disks from Fischer et al. (2016).

WISE_magnitude_distribution.png; dust_distance_V2.py; plots of distribution of magnitudes of each WISE band at a given distance

distance_v_temperature.png; dust_distance_V2.py; plot of distribution of temperatures and corresponding temperatures for each band (JHK & WISE)

~distribution_of_lifelines.png; uncertain_lifelines.py; plot of lifelines for samples generated through gaussian distribution around each object's ALMA-band flux and distance

~fitting_alpha_to_data.png; operation_Find_Alpha.py; plots of lifelines given false data and different α parameters under the model N(F) = A*F^-α with an overlay of the lifelines from real data.

light_curves.png; dust_distance_V2.py; plot of (normalized) light curves as flux vs. distance of disks for which there was detection across all bands (JHK & WISE). Temperature of binary stars are color-coded.

log_g_vs_Teff_post-cut.png; log_g_Teff_cuts.py; plot of log(g) vs. T-eff for binaries that  made the cut. The 'cut' being the log(g) vs. T-eff cutoff lines found by Kounkel et al. (2019) for cutting out objects not within the Orion cluster. These lines are also plotted.

log_g_vs_Teff_pre-cut.png; log_g_Teff_cuts.py; plot of log(g) vs. T-eff for all binaries with the Kounkel et al. (2019) cutoff lines also plotted. Error bars found from Kounkel et  al. (2019).

~vary_detection_mask.png; operation_FindAlpha.py; plots of lifelines for false data with randomized detection masks for one α value