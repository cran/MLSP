# MLSP 0.2.0

## New features

* `ml_f()` now compares every combination of preprocessing and calibration
  method in one call. By default it runs all 26 distinct preprocessing methods
  (raw reflectance, standard normal variate, absorbance, and Savitzky-Golay
  smoothing, first and second derivatives) with PCR, PLSR, LASSO, 'Cubist',
  and random forest, and returns the top 10 combinations. New arguments
  `models`, `top_n`, `rank_by`, and `verbose`; `smoother_selection` accepts
  `"all"`, method names, or a function.
* 26 preprocessing methods are available and compared by default: raw
  reflectance, SNV, log10(1/R), and 23 Savitzky-Golay filters (linear and
  quadratic smoothers, quadratic and cubic first derivatives, and quadratic
  second derivatives with 3- to 11-point windows).
* The 14 Savitzky-Golay settings of version 0.1.0 that reproduce another
  method exactly (cubic smoothers, linear first derivatives, cubic second
  derivatives, and the 3-point quadratic smoother) were removed; requesting
  one of them names the identical method to use instead.
* New `preprocess_methods()` lists the methods and explains how window size,
  polynomial order, and derivative order work. The methods and naming follow
  Clingensmith et al. (2019), doi:10.1111/ejss.12753.
* New `spectral_preprocess()` applies any method by name.
* `merge_of_lab_and_spectrum()` accepts any number of soil property columns
  and gains `soil_id`, `spec_id`, `n_splits`, and `cal_prop` arguments.
* `type_of_soil` accepts a column name as well as an index.

## Bug fixes

* Savitzky-Golay preprocessing failed with current versions of 'gsignal'
  (named spectra are no longer passed to `gsignal::sgolayfilt()`).
* Quadratic smoothing with 9- and 11-point windows returned the wrong object.
* PCR and PLSR predictions used all components at once, giving invalid
  metrics; the number of components is now chosen by cross-validation
  (one-standard-error rule) and reported in the `LV` column.
* Calibration metrics are now cross-validated for every model (out-of-bag for
  random forest); the `cal-mse` and `cal-rpiq` columns are renamed `cv-mse`
  and `cv-rpiq`.
* Random forest results were mislabelled by number of trees; `rf_ntree`
  (default 500) now sets the number of trees.
* `msd.comp()` no longer warns when predictions are constant.

## Other changes

* Maintainer changed to Katsutoshi Mizuta.
* All examples run on simulated spectra.
* Added a 'testthat' test suite covering every preprocessing method and the
  bug fixes above.
