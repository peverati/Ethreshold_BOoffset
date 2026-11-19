# Result: complete (J, v) threshold map of {X+, X+, e-}, M_c <= ~100 (threshold_map)

*Last updated: 2026-11-17 09:48*

Written by `analyze_map.py` from the `run_map.py` CSV files in `data/`. Crossings M* of fixed bases are rigorous upper bounds to M_c (a variational root below E_th proves binding), up to rounding. Values marked 'extrapolated' assume geometric convergence in the compact size.

## Blind-test integrity

* sha256(PREDICTION_MAP.md) = 6e494c54232cbadcc3423309babc48fbc8fd7a1a83a5cdb2f3fff8ef99e02fba: matches prediction_map.sha256.

## (2,1): L=2 'S', root 1

BO 29.019661, BO+ad 28.066887.

| basis | size | kept | M* (upper bound) | E_(k-1) - E_th at M* |
|---|---|---|---|---|
| compact N=12 | 1414 | 1291 | 27.7028828874 | -0.03220803812938322 |
| compact N=14 | 2104 | 1885 | 27.7000810601 | -0.03220433970279185 |
| compact N=16 | 2988 | 2630 | 27.6988979831 | -0.032202777894758006 |
| compact N=18 | 4090 | 3550 | 27.6982319738 | -0.032201898643480054 |
| compact N=16 + 1 diffuse | 3768 | 2942 | 27.6972832172 | -0.032200646069823646 |
| compact N=16 + 2 diffuse | 4788 | 3362 | 27.6972471656 | -0.03220059847237028 |
| compact N=16 + 3 diffuse | 6216 | 4061 | 27.6972470555 | -0.03220059832711847 |
| compact N=14 + 3 diffuse | 5332 | 3368 | 27.6972470559 | -0.03220059832746097 |
| compact N=18 + 3 diffuse | 7318 | 4930 | 27.6972470555 | -0.03220059832707528 |

* root ordering: roots 0..0 bound at every crossing (root k = v confirmed)
* **best proven upper bound (up to rounding): 27.697247055**
* rounding-safe bound (margin 1e-11 Eh): 27.697247072
* tol 1e-10 on 'compact N=16 + 3 diffuse': 27.6972470591 (difference from best: +3.6e-09)
* compact-size sequence (3 diffuse sets), N = [14, 16, 18]: 27.697247056, 27.697247056, 27.697247055; ratio 0.086
* **extrapolated M_c = 27.69724706 +/- 1.0e-09** (geometric in N; the bound is rigorous)
* diffuse increments at N_D: 0:27.6988979831, 1:27.6972832172, 2:27.6972471656, 3:27.6972470555; last increment -1.1e-07

| M | x = M/M* - 1 | E_th - E | E_b/(M - M*) |
|---|---|---|---|
| 27.71109568 | 0.0005 | 8.247017e-06 | 5.9551e-04 |
| 27.72494430 | 0.0010 | 1.652473e-05 | 5.9662e-04 |
| 27.75264155 | 0.0020 | 3.316406e-05 | 5.9869e-04 |
| 27.83573329 | 0.0050 | 8.367558e-05 | 6.0422e-04 |
| 27.97421953 | 0.0100 | 1.695506e-04 | 6.1216e-04 |
| 28.25119200 | 0.0200 | 3.464555e-04 | 6.2543e-04 |
| 29.08210941 | 0.0500 | 9.070311e-04 | 6.5496e-04 |
| 30.46697176 | 0.1000 | 1.904152e-03 | 6.8749e-04 |
| 33.23669647 | 0.2000 | 4.006045e-03 | 7.2318e-04 |

## (1,2): L=1 'A', root 2

BO 40.755428, BO+ad 39.766531.

| basis | size | kept | M* (upper bound) | E_(k-1) - E_th at M* |
|---|---|---|---|---|
| compact N=12 | 1365 | 1083 | 39.5256494848 | -0.018472961666730625 |
| compact N=14 | 2040 | 1569 | 39.4520098192 | -0.01842143767245147 |
| compact N=16 | 2907 | 2180 | 39.4109156775 | -0.018392647850849098 |
| compact N=18 | 3990 | 2919 | 39.3836265726 | -0.018373514856110607 |
| compact N=16 + 1 diffuse | 3297 | 2386 | 39.3191368182 | -0.018328252877631124 |
| compact N=16 + 2 diffuse | 3807 | 2667 | 39.3089487800 | -0.018321096408681348 |
| compact N=16 + 3 diffuse | 4521 | 3145 | 39.3086361518 | -0.018320876780708595 |
| compact N=16 + 4 diffuse | 5403 | 3585 | 39.3086244435 | -0.018320868555326775 |
| compact N=16 + 5 diffuse | 6579 | 4273 | 39.3086240179 | -0.018320868256350598 |
| compact N=14 + 5 diffuse | 5712 | 3709 | 39.3086240190 | -0.01832086825647744 |
| compact N=18 + 5 diffuse | 7662 | 4983 | 39.3086240178 | -0.018320868256386125 |

* root ordering: roots 0..1 bound at every crossing (root k = v confirmed)
* **best proven upper bound (up to rounding): 39.308624018**
* rounding-safe bound (margin 1e-11 Eh): 39.308624131
* tol 1e-10 on 'compact N=16 + 5 diffuse': 39.3086240281 (difference from best: +1.0e-08)
* compact-size sequence (5 diffuse sets), N = [14, 16, 18]: 39.308624019, 39.308624018, 39.308624018; ratio 0.098
* **extrapolated M_c = 39.30862402 +/- 1.0e-09** (geometric in N; the bound is rigorous)
* diffuse increments at N_D: 0:39.4109156775, 1:39.3191368182, 2:39.3089487800, 3:39.3086361518, 4:39.3086244435, 5:39.3086240179; last increment -4.3e-07

| M | x = M/M* - 1 | E_th - E | E_b/(M - M*) |
|---|---|---|---|
| 39.32827833 | 0.0005 | 1.855113e-06 | 9.4387e-05 |
| 39.34793264 | 0.0010 | 3.807602e-06 | 9.6864e-05 |
| 39.38724127 | 0.0020 | 7.898580e-06 | 1.0047e-04 |
| 39.50516714 | 0.0050 | 2.121791e-05 | 1.0796e-04 |
| 39.70171026 | 0.0100 | 4.594839e-05 | 1.1689e-04 |
| 40.09479650 | 0.0200 | 1.024621e-04 | 1.3033e-04 |
| 41.27405522 | 0.0500 | 3.126950e-04 | 1.5910e-04 |
| 43.23948642 | 0.1000 | 7.587020e-04 | 1.9301e-04 |
| 47.17034882 | 0.2000 | 1.876256e-03 | 2.3866e-04 |

## (3,1): L=3 'A', root 1

BO 42.355695, BO+ad 41.408480.

| basis | size | kept | M* (upper bound) | E_(k-1) - E_th at M* |
|---|---|---|---|---|
| compact N=10 | 1144 | 1002 | 41.0790339043 | -0.0269762246148032 |
| compact N=12 | 1820 | 1543 | 41.0763504962 | -0.026973542561598618 |
| compact N=14 | 2720 | 2243 | 41.0752485402 | -0.026972441086414833 |
| compact N=16 | 3876 | 3124 | 41.0748086490 | -0.0269720013773756 |
| compact N=14 + 1 diffuse | 3500 | 2776 | 41.0743475121 | -0.026971540423364748 |
| compact N=14 + 2 diffuse | 4520 | 3345 | 41.0743472743 | -0.026971540185309284 |
| compact N=14 + 3 diffuse | 5948 | 4279 | 41.0743472675 | -0.026971540178698072 |
| compact N=12 + 3 diffuse | 5048 | 3629 | 41.0743472808 | -0.02697154019090925 |
| compact N=16 + 3 diffuse | 7104 | 5110 | 41.0743472647 | -0.02697154017678549 |

* root ordering: roots 0..0 bound at every crossing (root k = v confirmed)
* **best proven upper bound (up to rounding): 41.074347265**
* rounding-safe bound (margin 1e-11 Eh): 41.074347285
* tol 1e-10 on 'compact N=14 + 3 diffuse': 41.0743477125 (difference from best: +4.5e-07)
* compact-size sequence (3 diffuse sets), N = [12, 14, 16]: 41.074347281, 41.074347268, 41.074347265; ratio 0.215
* **extrapolated M_c = 41.07434726 +/- 1.0e-09** (geometric in N; the bound is rigorous)
* diffuse increments at N_D: 0:41.0752485402, 1:41.0743475121, 2:41.0743472743, 3:41.0743472675; last increment -6.8e-09

| M | x = M/M* - 1 | E_th - E | E_b/(M - M*) |
|---|---|---|---|
| 41.09488444 | 0.0005 | 1.173121e-05 | 5.7122e-04 |
| 41.11542161 | 0.0010 | 2.347669e-05 | 5.7157e-04 |
| 41.15649596 | 0.0020 | 4.700934e-05 | 5.7225e-04 |
| 41.27971900 | 0.0050 | 1.179242e-04 | 5.7420e-04 |
| 41.48509074 | 0.0100 | 2.370832e-04 | 5.7720e-04 |
| 41.89583421 | 0.0200 | 4.785321e-04 | 5.8252e-04 |
| 43.12806463 | 0.0500 | 1.221462e-03 | 5.9476e-04 |
| 45.18178199 | 0.1000 | 2.495886e-03 | 6.0765e-04 |
| 49.28921672 | 0.2000 | 5.079565e-03 | 6.1834e-04 |

## (5,0): L=5 'A', root 0

BO 50.318074, BO+ad 49.353402.

| basis | size | kept | M* (upper bound) | E_(k-1) - E_th at M* |
|---|---|---|---|---|
| compact N=8 | 990 | 893 | 49.1154786408 | nan |
| compact N=10 | 1716 | 1473 | 49.1153834277 | nan |
| compact N=12 | 2730 | 2238 | 49.1153782330 | nan |
| compact N=14 | 4080 | 3214 | 49.1153751438 | nan |
| compact N=12 + 1 diffuse | 3900 | 3017 | 49.1153808129 | nan |
| compact N=12 + 2 diffuse | 5430 | 3855 | 49.1153803378 | nan |
| compact N=10 + 2 diffuse | 4416 | 3177 | 49.1153830092 | nan |
| compact N=14 + 2 diffuse | 6780 | 4744 | 49.1153754484 | nan |
* **best proven upper bound (up to rounding): 49.115375144**
* rounding-safe bound (margin 1e-11 Eh): 49.115380349
* tol 1e-10 on 'compact N=12 + 2 diffuse': 49.1155087239 (difference from best: +1.3e-04)
* compact-size sequence (2 diffuse sets), N = [10, 12, 14]: 49.115383009, 49.115380338, 49.115375448; ratio nan
* **extrapolated M_c = 49.11537545 +/- 4.9e-06** (geometric in N; the bound is rigorous)
* diffuse increments at N_D: 0:49.1153782330, 1:49.1153808129, 2:49.1153803378; last increment -4.8e-07

| M | x = M/M* - 1 | E_th - E | E_b/(M - M*) |
|---|---|---|---|
| 49.13993283 | 0.0005 | 2.196335e-05 | 8.9436e-04 |
| 49.16449052 | 0.0010 | 4.392767e-05 | 8.9438e-04 |
| 49.21360589 | 0.0020 | 8.784515e-05 | 8.9427e-04 |
| 49.36095202 | 0.0050 | 2.195059e-04 | 8.9384e-04 |
| 49.60652890 | 0.0100 | 4.386199e-04 | 8.9304e-04 |
| 50.09768265 | 0.0200 | 8.755453e-04 | 8.9131e-04 |
| 51.57114390 | 0.0500 | 2.174280e-03 | 8.8538e-04 |
| 54.02691266 | 0.1000 | 4.291355e-03 | 8.7373e-04 |
| 58.93845017 | 0.2000 | 8.320495e-03 | 8.4704e-04 |

## (2,2): L=2 'S', root 2

BO 53.968672, BO+ad 53.002721.

| basis | size | kept | M* (upper bound) | E_(k-1) - E_th at M* |
|---|---|---|---|---|
| compact N=10 | 1341 | 984 | 52.7335617779 | -0.01761638681742661 |
| compact N=12 | 2121 | 1492 | 52.6668204209 | -0.017576618663166843 |
| compact N=14 | 3156 | 2146 | 52.6375312765 | -0.017559151734060918 |
| compact N=16 | 4482 | 2959 | 52.6173367046 | -0.017547103185889756 |
| compact N=14 + 1 diffuse | 3936 | 2543 | 52.5892741745 | -0.017530353281352384 |
| compact N=14 + 2 diffuse | 4956 | 2970 | 52.5891648488 | -0.01753028801098655 |
| compact N=14 + 3 diffuse | 6384 | 3685 | 52.5891644775 | -0.017530287789297494 |
| compact N=12 + 3 diffuse | 5349 | 3063 | 52.5891645291 | -0.01753028781798366 |
| compact N=16 + 3 diffuse | 7710 | 4470 | 52.5891644679 | -0.017530287783766196 |

* root ordering: roots 0..1 bound at every crossing (root k = v confirmed)
* **best proven upper bound (up to rounding): 52.589164468**
* rounding-safe bound (margin 1e-11 Eh): 52.589164528
* tol 1e-10 on 'compact N=14 + 3 diffuse': 52.5891645171 (difference from best: +4.9e-08)
* compact-size sequence (3 diffuse sets), N = [12, 14, 16]: 52.589164529, 52.589164477, 52.589164468; ratio 0.185
* **extrapolated M_c = 52.58916447 +/- 2.2e-09** (geometric in N; the bound is rigorous)
* diffuse increments at N_D: 0:52.6375312765, 1:52.5892741745, 2:52.5891648488, 3:52.5891644775; last increment -3.7e-07

| M | x = M/M* - 1 | E_th - E | E_b/(M - M*) |
|---|---|---|---|
| 52.61545905 | 0.0005 | 5.192028e-06 | 1.9746e-04 |
| 52.64175363 | 0.0010 | 1.041846e-05 | 1.9811e-04 |
| 52.69434280 | 0.0020 | 2.096477e-05 | 1.9933e-04 |
| 52.85211029 | 0.0050 | 5.326269e-05 | 2.0256e-04 |
| 53.11505611 | 0.0100 | 1.089738e-04 | 2.0722e-04 |
| 53.64094776 | 0.0200 | 2.262206e-04 | 2.1508e-04 |
| 55.21862269 | 0.0500 | 6.131233e-04 | 2.3317e-04 |
| 57.84808091 | 0.1000 | 1.338571e-03 | 2.5453e-04 |
| 63.10699736 | 0.2000 | 2.959616e-03 | 2.8139e-04 |

## (4,1): L=4 'S', root 1

BO 58.421973, BO+ad 57.475700.

| basis | size | kept | M* (upper bound) | E_(k-1) - E_th at M* |
|---|---|---|---|---|
| compact N=10 | 1466 | 1219 | 57.1665921360 | -0.02309191072503214 |
| compact N=12 | 2324 | 1853 | 57.1655029775 | -0.023091072241361366 |
| compact N=14 | 3464 | 2675 | 57.1651674632 | -0.023090813842412128 |
| compact N=16 | 4926 | 3694 | 57.1650424902 | -0.023090717577921194 |
| compact N=14 + 1 diffuse | 4634 | 3388 | 57.1649492291 | -0.023090645789004827 |
| compact N=14 + 2 diffuse | 6164 | 4098 | 57.1649492270 | -0.023090645788119424 |
| compact N=12 + 2 diffuse | 5024 | 3328 | 57.1649490696 | -0.0230906457595747 |
| compact N=16 + 2 diffuse | 7626 | 5067 | 57.1649493512 | -0.023090645814119626 |

* root ordering: roots 0..0 bound at every crossing (root k = v confirmed)
* **best proven upper bound (up to rounding): 57.164949070**
* rounding-safe bound (margin 1e-11 Eh): 57.164949247
* tol 1e-10 on 'compact N=14 + 2 diffuse': 57.1649984671 (difference from best: +4.9e-05)
* compact-size sequence (2 diffuse sets), N = [12, 14, 16]: 57.164949070, 57.164949227, 57.164949351; ratio 0.789
* **extrapolated M_c = 57.16494981 +/- 4.6e-07** (geometric in N; the bound is rigorous)
* diffuse increments at N_D: 0:57.1651674632, 1:57.1649492291, 2:57.1649492270; last increment -2.0e-09

| M | x = M/M* - 1 | E_th - E | E_b/(M - M*) |
|---|---|---|---|
| 57.19353154 | 0.0005 | 1.421402e-05 | 4.9730e-04 |
| 57.22211402 | 0.0010 | 2.843631e-05 | 4.9744e-04 |
| 57.27927897 | 0.0020 | 5.690503e-05 | 4.9773e-04 |
| 57.45077381 | 0.0050 | 1.424968e-04 | 4.9855e-04 |
| 57.73659856 | 0.0100 | 2.857234e-04 | 4.9982e-04 |
| 58.30824805 | 0.0200 | 5.740388e-04 | 5.0209e-04 |
| 60.02319652 | 0.0500 | 1.449535e-03 | 5.0714e-04 |
| 62.88144398 | 0.1000 | 2.924517e-03 | 5.1159e-04 |
| 68.59793888 | 0.2000 | 5.852967e-03 | 5.1194e-04 |

## (6,0): L=6 'S', root 0

BO 68.181554, BO+ad 67.215424.

| basis | size | kept | M* (upper bound) | E_(k-1) - E_th at M* |
|---|---|---|---|---|
| compact N=6 | 604 | 545 | 66.9843516080 | nan |
| compact N=8 | 1180 | 997 | 66.9838176233 | nan |
| compact N=10 | 2038 | 1618 | 66.9837924238 | nan |
| compact N=12 | 3234 | 2433 | 66.9837577570 | nan |
| compact N=10 + 1 diffuse | 3598 | 2637 | 66.9838010195 | nan |
| compact N=10 + 2 diffuse | 5638 | 3617 | 66.9837986187 | nan |
| compact N=8 + 2 diffuse | 4780 | 3074 | 66.9838268481 | nan |
| compact N=12 + 2 diffuse | 6834 | 4353 | 66.9837745786 | nan |
* **best proven upper bound (up to rounding): 66.983757757**
* rounding-safe bound (margin 1e-11 Eh): 66.983798633
* tol 1e-10 on 'compact N=10 + 2 diffuse': 66.9882392562 (difference from best: +4.5e-03)
* compact-size sequence (2 diffuse sets), N = [8, 10, 12]: 66.983826848, 66.983798619, 66.983774579; ratio 0.852
* **extrapolated M_c = 66.98363663 +/- 1.4e-04** (geometric in N; the bound is rigorous)
* diffuse increments at N_D: 0:66.9837924238, 1:66.9838010195, 2:66.9837986187; last increment -2.4e-06

| M | x = M/M* - 1 | E_th - E | E_b/(M - M*) |
|---|---|---|---|
| 67.01724964 | 0.0005 | 2.291204e-05 | 6.8411e-04 |
| 67.05074151 | 0.0010 | 4.584671e-05 | 6.8445e-04 |
| 67.11772527 | 0.0020 | 9.169987e-05 | 6.8449e-04 |
| 67.31867655 | 0.0050 | 2.291282e-04 | 6.8413e-04 |
| 67.65359533 | 0.0100 | 4.577275e-04 | 6.8334e-04 |
| 68.32343291 | 0.0200 | 9.131714e-04 | 6.8164e-04 |
| 70.33294564 | 0.0500 | 2.264249e-03 | 6.7606e-04 |
| 73.68213353 | 0.1000 | 4.459683e-03 | 6.6579e-04 |
| 80.38050931 | 0.2000 | 8.620328e-03 | 6.4346e-04 |

## (3,2): L=3 'A', root 2

BO 70.347388, BO+ad 69.393851.

| basis | size | kept | M* (upper bound) | E_(k-1) - E_th at M* |
|---|---|---|---|---|
| compact N=8 | 990 | 751 | 69.1008604626 | -0.016332322154652135 |
| compact N=10 | 1716 | 1222 | 69.0456598649 | -0.016304720061115052 |
| compact N=12 | 2730 | 1858 | 69.0278679645 | -0.01629581831911897 |
| compact N=14 | 4080 | 2662 | 69.0208233685 | -0.01629229309809671 |
| compact N=12 + 1 diffuse | 3510 | 2443 | 69.0138900372 | -0.016288823215829085 |
| compact N=12 + 2 diffuse | 4530 | 3017 | 69.0138890647 | -0.01628882272939841 |
| compact N=12 + 3 diffuse | 5958 | 3954 | 69.0138887790 | -0.016288822586582208 |
| compact N=10 + 3 diffuse | 4944 | 3364 | 69.0138895184 | -0.016288822883121612 |
| compact N=14 + 3 diffuse | 7308 | 4721 | 69.0138886605 | -0.01628882253112629 |

* root ordering: roots 0..1 bound at every crossing (root k = v confirmed)
* **best proven upper bound (up to rounding): 69.013888660**
* rounding-safe bound (margin 1e-11 Eh): 69.013888822
* tol 1e-10 on 'compact N=12 + 3 diffuse': 69.0138902395 (difference from best: +1.6e-06)
* compact-size sequence (3 diffuse sets), N = [10, 12, 14]: 69.013889518, 69.013888779, 69.013888660; ratio 0.160
* **extrapolated M_c = 69.01388864 +/- 2.3e-08** (geometric in N; the bound is rigorous)
* diffuse increments at N_D: 0:69.0278679645, 1:69.0138900372, 2:69.0138890647, 3:69.0138887790; last increment -2.9e-07

| M | x = M/M* - 1 | E_th - E | E_b/(M - M*) |
|---|---|---|---|
| 69.04839560 | 0.0005 | 8.021337e-06 | 2.3246e-04 |
| 69.08290255 | 0.0010 | 1.606210e-05 | 2.3274e-04 |
| 69.15191644 | 0.0020 | 3.220028e-05 | 2.3329e-04 |
| 69.35895810 | 0.0050 | 8.104564e-05 | 2.3487e-04 |
| 69.70402755 | 0.0100 | 1.637772e-04 | 2.3731e-04 |
| 70.39416643 | 0.0200 | 3.335934e-04 | 2.4169e-04 |
| 72.46458309 | 0.0500 | 8.703351e-04 | 2.5222e-04 |
| 75.91527753 | 0.1000 | 1.826594e-03 | 2.6467e-04 |
| 82.81666639 | 0.2000 | 3.853776e-03 | 2.7920e-04 |

## (1,3): L=1 'A', root 3

BO 72.664208, BO+ad 71.661196.

| basis | size | kept | M* (upper bound) | E_(k-1) - E_th at M* |
|---|---|---|---|---|
| compact N=12 | 1365 | 1050 | 72.2149350093 | -0.010582093805915849 |
| compact N=14 | 2040 | 1520 | 71.9183894750 | -0.010481849092746354 |
| compact N=16 | 2907 | 2107 | 71.7315877359 | -0.01041859449870114 |
| compact N=18 | 3990 | 2827 | 71.6006648249 | -0.01037421213543882 |
| compact N=16 + 1 diffuse | 3297 | 2353 | 71.2132283560 | -0.010242635829532387 |
| compact N=16 + 2 diffuse | 3807 | 2639 | 71.1918719967 | -0.010235372811997434 |
| compact N=16 + 3 diffuse | 4521 | 3120 | 71.1912559053 | -0.010235163271568037 |
| compact N=16 + 4 diffuse | 5403 | 3567 | 71.1912330523 | -0.010235155498944104 |
| compact N=16 + 5 diffuse | 6579 | 4258 | 71.1912322178 | -0.010235155215169989 |
| compact N=14 + 5 diffuse | 5712 | 3691 | 71.1912322302 | -0.010235155218309921 |
| compact N=18 + 5 diffuse | 7662 | 4953 | 71.1912322167 | -0.010235155214927405 |

* root ordering: roots 0..2 bound at every crossing (root k = v confirmed)
* **best proven upper bound (up to rounding): 71.191232217**
* rounding-safe bound (margin 1e-11 Eh): 71.191232553
* tol 1e-10 on 'compact N=16 + 5 diffuse': 71.1912322365 (difference from best: +2.0e-08)
* compact-size sequence (5 diffuse sets), N = [14, 16, 18]: 71.191232230, 71.191232218, 71.191232217; ratio 0.090
* **extrapolated M_c = 71.19123222 +/- 1.0e-09** (geometric in N; the bound is rigorous)
* diffuse increments at N_D: 0:71.7315877359, 1:71.2132283560, 2:71.1918719967, 3:71.1912559053, 4:71.1912330523, 5:71.1912322178; last increment -8.3e-07

| M | x = M/M* - 1 | E_th - E | E_b/(M - M*) |
|---|---|---|---|
| 71.22682783 | 0.0005 | 1.146575e-06 | 3.2211e-05 |
| 71.26242345 | 0.0010 | 2.369107e-06 | 3.3278e-05 |
| 71.33361468 | 0.0020 | 4.962441e-06 | 3.4853e-05 |
| 71.54718838 | 0.0050 | 1.359844e-05 | 3.8203e-05 |
| 71.90314454 | 0.0100 | 3.013391e-05 | 4.2328e-05 |
| 72.61505686 | 0.0200 | 6.943611e-05 | 4.8767e-05 |
| 74.75079383 | 0.0500 | 2.255275e-04 | 6.3358e-05 |
| 78.31035544 | 0.1000 | 5.816652e-04 | 8.1705e-05 |
| 85.42947866 | 0.2000 | 1.537144e-03 | 1.0796e-04 |

## (5,1): L=5 'A', root 1

BO 77.189949, BO+ad 76.242851.

| basis | size | kept | M* (upper bound) | E_(k-1) - E_th at M* |
|---|---|---|---|---|
| compact N=8 | 990 | 825 | 75.9513502603 | -0.020142837565587213 |
| compact N=10 | 1716 | 1354 | 75.9492934850 | -0.02014159570657359 |
| compact N=12 | 2730 | 2046 | 75.9488413631 | -0.020141321478548357 |
| compact N=14 | 4080 | 2934 | 75.9487150938 | -0.02014124499952963 |
| compact N=12 + 1 diffuse | 3900 | 2974 | 75.9486537585 | -0.02014120772359107 |
| compact N=12 + 2 diffuse | 5430 | 3824 | 75.9486539770 | -0.020141207767030378 |
| compact N=10 + 2 diffuse | 4416 | 3174 | 75.9486542961 | -0.020141208186305604 |
| compact N=14 + 2 diffuse | 6780 | 4670 | 75.9486537358 | -0.020141207623033175 |

* root ordering: roots 0..0 bound at every crossing (root k = v confirmed)
* **best proven upper bound (up to rounding): 75.948653736**
* rounding-safe bound (margin 1e-11 Eh): 75.948654001
* tol 1e-10 on 'compact N=12 + 2 diffuse': 75.9489625435 (difference from best: +3.1e-04)
* compact-size sequence (2 diffuse sets), N = [10, 12, 14]: 75.948654296, 75.948653977, 75.948653736; ratio 0.756
* **extrapolated M_c = 75.94865299 +/- 7.5e-07** (geometric in N; the bound is rigorous)
* diffuse increments at N_D: 0:75.9488413631, 1:75.9486537585, 2:75.9486539770; last increment +2.2e-07

| M | x = M/M* - 1 | E_th - E | E_b/(M - M*) |
|---|---|---|---|
| 75.98662806 | 0.0005 | 1.605469e-05 | 4.2278e-04 |
| 76.02460239 | 0.0010 | 3.211418e-05 | 4.2284e-04 |
| 76.10055104 | 0.0020 | 6.424696e-05 | 4.2296e-04 |
| 76.32839700 | 0.0050 | 1.607510e-04 | 4.2331e-04 |
| 76.70814027 | 0.0100 | 3.219127e-04 | 4.2386e-04 |
| 77.46762681 | 0.0200 | 6.452386e-04 | 4.2479e-04 |
| 79.74608642 | 0.0500 | 1.619954e-03 | 4.2659e-04 |
| 83.54351911 | 0.1000 | 3.244850e-03 | 4.2724e-04 |
| 91.13838448 | 0.2000 | 6.429730e-03 | 4.2329e-04 |

## (2,3): L=2 'S', root 3

BO 88.206149, BO+ad 87.225441.

| basis | size | kept | M* (upper bound) | E_(k-1) - E_th at M* |
|---|---|---|---|---|
| compact N=10 | 1341 | 1007 | 87.4054857701 | -0.010857972246653469 |
| compact N=12 | 2121 | 1531 | 87.1485745571 | -0.010778531454334717 |
| compact N=14 | 3156 | 2204 | 86.9882221620 | -0.010728877972301099 |
| compact N=16 | 4482 | 3040 | 86.9135852326 | -0.010705748251708003 |
| compact N=14 + 1 diffuse | 3936 | 2629 | 86.7875790541 | -0.01066667317554959 |
| compact N=14 + 2 diffuse | 4956 | 3061 | 86.7873322225 | -0.010666596599680278 |
| compact N=14 + 3 diffuse | 6384 | 3776 | 86.7873306218 | -0.010666596103101889 |
| compact N=12 + 3 diffuse | 5349 | 3135 | 86.7873310457 | -0.010666596215016755 |
| compact N=16 + 3 diffuse | 7710 | 4582 | 86.7873305255 | -0.010666596074176637 |

* root ordering: roots 0..2 bound at every crossing (root k = v confirmed)
* **best proven upper bound (up to rounding): 86.787330525**
* rounding-safe bound (margin 1e-11 Eh): 86.787330745
* tol 1e-10 on 'compact N=14 + 3 diffuse': 86.7873308317 (difference from best: +3.1e-07)
* compact-size sequence (3 diffuse sets), N = [12, 14, 16]: 86.787331046, 86.787330622, 86.787330525; ratio 0.227
* **extrapolated M_c = 86.78733050 +/- 2.8e-08** (geometric in N; the bound is rigorous)
* diffuse increments at N_D: 0:86.9882221620, 1:86.7875790541, 2:86.7873322225, 3:86.7873306218; last increment -1.6e-06

| M | x = M/M* - 1 | E_th - E | E_b/(M - M*) |
|---|---|---|---|
| 86.83072419 | 0.0005 | 3.531721e-06 | 8.1388e-05 |
| 86.87411786 | 0.0010 | 7.098363e-06 | 8.1790e-05 |
| 86.96090519 | 0.0020 | 1.432590e-05 | 8.2535e-05 |
| 87.22126718 | 0.0050 | 3.667129e-05 | 8.4508e-05 |
| 87.65520383 | 0.0100 | 7.580901e-05 | 8.7350e-05 |
| 88.52307714 | 0.0200 | 1.600117e-04 | 9.2186e-05 |
| 91.12669705 | 0.0500 | 4.493735e-04 | 1.0356e-04 |
| 95.46606358 | 0.1000 | 1.019997e-03 | 1.1753e-04 |
| 104.14479663 | 0.2000 | 2.364486e-03 | 1.3622e-04 |

## (7,0): L=7 'A', root 0

BO 88.715954, BO+ad 87.748556.

**Computed with double-double matrix assembly** (paper1/threshold_map_L7hp, run_map_hp.py; validate_L 7 PASSES there, validate_hp.txt). Standard projection tol 1e-14, check tol 1e-12; otherwise run 5's method. (Run 5's double-precision path failed validate_L 7.)

Excluded rows (spurious-root guard / conditioning, recheck_guard.py):

| basis | tol | margin | M* | reason |
|---|---|---|---|---|
| compact N=8 + 2 diffuse | 1e-14 | 0.0 | nan | bound at the lower bracket end M=86.2486 (gap -1.535e-05): spurious level; conditioning: lam_min/lam_max=-3.1e-16 < -1e-16 |

| basis | size | kept | M* (upper bound) | E_(k-1) - E_th at M* |
|---|---|---|---|---|
| compact N=6 | 672 | 644 | 87.5224237646 | nan |
| compact N=8 | 1320 | 1191 | 87.5215195891 | nan |
| compact N=10 | 2288 | 1955 | 87.5214750528 | nan |
| compact N=12 | 3640 | 2953 | 87.5214699213 | nan |
| compact N=10 + 1 diffuse | 3848 | 3269 | 87.5214716945 | nan |
| compact N=10 + 2 diffuse | 5888 | 4481 | 87.5214715905 | nan |
| compact N=8 + 2 diffuse | 4920 | 3785 | nan | nan |
| compact N=12 + 2 diffuse | 7240 | 5409 | 87.5214689350 | nan |
* **best proven upper bound (up to rounding): 87.521468935**
* rounding-safe bound (margin 1e-11 Eh): 87.521471609
* tol 1e-12 on 'compact N=10 + 2 diffuse': 87.5218653241 (difference from best: +4.0e-04)
* diffuse increments at N_D: 0:87.5214750528, 1:87.5214716945, 2:87.5214715905; last increment -1.0e-07

| M | x = M/M* - 1 | E_th - E | E_b/(M - M*) |
|---|---|---|---|
| 87.56522967 | 0.0005 | 2.364274e-05 | 5.4027e-04 |
| 87.60899040 | 0.0010 | 4.728040e-05 | 5.4021e-04 |
| 87.69651187 | 0.0020 | 9.453606e-05 | 5.4007e-04 |
| 87.95907628 | 0.0050 | 2.361445e-04 | 5.3963e-04 |
| 88.39668362 | 0.0100 | 4.716216e-04 | 5.3886e-04 |
| 89.27189831 | 0.0200 | 9.405029e-04 | 5.3730e-04 |
| 91.89754238 | 0.0500 | 2.329604e-03 | 5.3235e-04 |
| 96.27361583 | 0.1000 | 4.582004e-03 | 5.2353e-04 |
| 105.02576272 | 0.2000 | 8.838273e-03 | 5.0492e-04 |

J7v0 summary (threshold_map_L7hp, 2026-11-17):

* Accepted crossings (tol 1e-14, double-double assembly, single-pass orthogonalization; every row passed run 6's crossing-slope guard (c): gap(M*(1-1e-3)) ~ +1.5e-5..+4.7e-5, gap(M*(1+1e-3)) = -4.7e-5): compact N = 6..12: 87.5224238, 87.5215196, 87.5214751, 87.5214699 (monotone); N = 10 + 1/2 diffuse: 87.5214717, 87.5214716; N = 12 + 2 diffuse: 87.5214689.
* Excluded: compact N = 8 + 2 diffuse (tol 1e-14): spurious level (bound at M = 86.25 although it is a subset of N = 10 + 2 diffuse; lam_min/lam_max = -3.1e-16). Cause: single-pass canonical orthogonalization at tol 1e-14 (diag_spurious.txt). With a second orthogonalization pass the same basis gives a normal crossing, 87.5215171.
* Numerical cross-check: a complete rerun with two-pass orthogonalization (data_2pass/) reproduces every accepted crossing to <= 1e-7 for N >= 10 (3e-6 at N = 6). Two-pass is not used for the results (it fails validate_hp at tol 1e-14/1e-15).
* **M_c(J7v0) = 87.521468 +/- 1e-6** (best bound 87.5214689 minus the geometric compact-size correction 6.7e-7, ratio 0.115; error = compact correction + last diffuse increment 1e-7 + numerical cross-check 1e-7). **Certified bound** (margin 1e-11 Eh at the standard tol 1e-14, N = 10 + 2 diffuse): 87.5214716. tol 1e-12 check: 87.5218653 (+3.9e-4, truncation; looser but valid bound).
* Blind test: frozen prediction 87.5223 +/- 0.0295 (PREDICTION_MAP.md, hash OK): exact - pred = -0.0008: **PASS**.

## (4,2): L=4 'S', root 2

BO 89.656384, BO+ad 88.708986.

| basis | size | kept | M* (upper bound) | E_(k-1) - E_th at M* |
|---|---|---|---|---|
| compact N=8 | 1275 | 893 | 88.3969287438 | -0.015033578674182368 |
| compact N=10 | 2199 | 1442 | 88.3681821888 | -0.015021534330523778 |
| compact N=12 | 3486 | 2169 | 88.3590716536 | -0.015017714776000168 |
| compact N=14 | 5196 | 3091 | 88.3561495103 | -0.015016489289837065 |
| compact N=12 + 1 diffuse | 4656 | 2928 | 88.3541940148 | -0.01501566952279193 |
| compact N=12 + 2 diffuse | 6186 | 3644 | 88.3541938075 | -0.015015669431745315 |
| compact N=10 + 2 diffuse | 4899 | 2960 | 88.3541954719 | -0.015015670107099255 |
| compact N=14 + 2 diffuse | 7896 | 4527 | 88.3541941296 | -0.015015669351430005 |

* root ordering: roots 0..1 bound at every crossing (root k = v confirmed)
* **best proven upper bound (up to rounding): 88.354193808**
* rounding-safe bound (margin 1e-11 Eh): 88.354193850
* tol 1e-10 on 'compact N=12 + 2 diffuse': 88.3543274154 (difference from best: +1.3e-04)
* compact-size sequence (2 diffuse sets), N = [10, 12, 14]: 88.354195472, 88.354193808, 88.354194130; ratio nan
* **extrapolated M_c = 88.35419413 +/- 3.2e-07** (geometric in N; the bound is rigorous)
* diffuse increments at N_D: 0:88.3590716536, 1:88.3541940148, 2:88.3541938075; last increment -2.1e-07

| M | x = M/M* - 1 | E_th - E | E_b/(M - M*) |
|---|---|---|---|
| 88.39837090 | 0.0005 | 1.029885e-05 | 2.3313e-04 |
| 88.44254800 | 0.0010 | 2.061120e-05 | 2.3328e-04 |
| 88.53090220 | 0.0020 | 4.127577e-05 | 2.3358e-04 |
| 88.79596478 | 0.0050 | 1.035782e-04 | 2.3446e-04 |
| 89.23773575 | 0.0100 | 2.083842e-04 | 2.3585e-04 |
| 90.12127768 | 0.0200 | 4.212524e-04 | 2.3839e-04 |
| 92.77190350 | 0.0500 | 1.080438e-03 | 2.4457e-04 |
| 97.18961319 | 0.1000 | 2.223772e-03 | 2.5169e-04 |
| 106.02503257 | 0.2000 | 4.576476e-03 | 2.5898e-04 |

## (6,1): L=6 'S', root 1

BO 98.647176, BO+ad 97.698595.

| basis | size | kept | M* (upper bound) | E_(k-1) - E_th at M* |
|---|---|---|---|---|
| compact N=6 | 604 | 473 | 97.4580753540 | -0.017858238049719255 |
| compact N=8 | 1180 | 852 | 97.4226508222 | -0.01784196810439581 |
| compact N=10 | 2038 | 1382 | 97.4180129377 | -0.017839761002402088 |
| compact N=12 | 3234 | 2078 | 97.4170898580 | -0.017839307171320262 |
| compact N=10 + 1 diffuse | 3598 | 2537 | 97.4167740427 | -0.017839157722751797 |
| compact N=10 + 2 diffuse | 5638 | 3531 | 97.4167766613 | -0.01783915870493591 |
| compact N=8 + 2 diffuse | 4780 | 3039 | 97.4172322483 | -0.017839328561497447 |
| compact N=12 + 2 diffuse | 6834 | 4186 | 97.4168010825 | -0.017839161092406408 |

* root ordering: roots 0..0 bound at every crossing (root k = v confirmed)
* **best proven upper bound (up to rounding): 97.416774043**
* rounding-safe bound (margin 1e-11 Eh): 97.416776689
* tol 1e-10 on 'compact N=10 + 2 diffuse': 97.4274183377 (difference from best: +1.1e-02)
* compact-size sequence (2 diffuse sets), N = [8, 10, 12]: 97.417232248, 97.416776661, 97.416801083; ratio nan
* **extrapolated M_c = 97.41680108 +/- 2.4e-05** (geometric in N; the bound is rigorous)
* diffuse increments at N_D: 0:97.4180129377, 1:97.4167740427, 2:97.4167766613; last increment +2.6e-06

| M | x = M/M* - 1 | E_th - E | E_b/(M - M*) |
|---|---|---|---|
| 97.46548243 | 0.0005 | 1.746701e-05 | 3.5860e-04 |
| 97.51419082 | 0.0010 | 3.493725e-05 | 3.5864e-04 |
| 97.61160759 | 0.0020 | 6.988444e-05 | 3.5869e-04 |
| 97.90385791 | 0.0050 | 1.747759e-04 | 3.5882e-04 |
| 98.39094178 | 0.0100 | 3.497382e-04 | 3.5901e-04 |
| 99.36510952 | 0.0200 | 7.000391e-04 | 3.5930e-04 |
| 102.28761274 | 0.0500 | 1.751302e-03 | 3.5955e-04 |
| 107.15845145 | 0.1000 | 3.491812e-03 | 3.5844e-04 |
| 116.90012885 | 0.2000 | 6.873883e-03 | 3.5281e-04 |

## Summary: new thresholds

| (J,v) | BO | BO+ad | best upper bound | extrapolated M_c | prediction (frozen) | exact - pred | pass | Frolov y_t (bound at) | consistent with Frolov bound |
|---|---|---|---|---|---|---|---|---|---|
| (2,1) | 29.01966 | 28.06689 | 27.6972471 | 27.6972471 | 27.7068 +/- 0.0295 | -0.0095 | PASS | 27.715 (27.75) | yes |
| (1,2) | 40.75543 | 39.76653 | 39.3086240 | 39.3086240 | 39.3093 +/- 0.0295 | -0.0007 | PASS | 39.375 (39.45) | yes |
| (3,1) | 42.35569 | 41.40848 | 41.0743473 | 41.0743473 | 41.0810 +/- 0.0295 | -0.0067 | PASS | | |
| (5,0) | 50.31807 | 49.35340 | 49.1153751 | 49.1153754 | 49.1173 +/- 0.0295 | -0.0019 | PASS | | |
| (2,2) | 53.96867 | 53.00272 | 52.5891645 | 52.5891645 | 52.5920 +/- 0.0295 | -0.0029 | PASS | 53.295 (53.3) | yes |
| (4,1) | 58.42197 | 57.47570 | 57.1649491 | 57.1649498 | 57.1697 +/- 0.0295 | -0.0047 | PASS | | |
| (6,0) | 68.18155 | 67.21542 | 66.9837578 | 66.9836366 | 66.9849 +/- 0.0295 | -0.0013 | PASS | | |
| (3,2) | 70.34739 | 69.39385 | 69.0138887 | 69.0138886 | 69.0171 +/- 0.0295 | -0.0032 | PASS | | |
| (1,3) | 72.66421 | 71.66120 | 71.1912322 | 71.1912322 | 71.1896 +/- 0.0295 | +0.0016 | PASS | | |
| (5,1) | 77.18995 | 76.24285 | 75.9486537 | 75.9486530 | 75.9521 +/- 0.0295 | -0.0034 | PASS | | |
| (2,3) | 88.20615 | 87.22544 | 86.7873305 | 86.7873305 | 86.7875 +/- 0.0295 | -0.0002 | PASS | | |
| (7,0) | 88.71595 | 87.74856 | 87.5214689 | 87.5214689 | 87.5223 +/- 0.0295 | -0.0008 | PASS | | |
| (4,2) | 89.65638 | 88.70899 | 88.3541938 | 88.3541941 | 88.3571 +/- 0.0295 | -0.0029 | PASS | | |
| (6,1) | 98.64718 | 97.69860 | 97.4167740 | 97.4168011 | 97.4193 +/- 0.0295 | -0.0025 | PASS | | |

Blind test: 14 PASS out of 14 targets with an exact value (criterion fixed in PREDICTION_MAP.md: |exact - pred| <= band + exact error).

## Complete map (M_c <= ~101)

| order | (J,v) | M_c | source | BO | BO+ad | BO - exact | BO+ad - exact |
|---|---|---|---|---|---|---|---|
| 1 | (1,0) | 4.208283 | exact (earlier) | 5.5375 | 4.5706 | 1.3292 | 0.3623 |
| 2 | (0,1) | 9.525634 | exact (earlier) | 11.0815 | 10.0632 | 1.5558 | 0.5375 |
| 3 | (2,0) | 11.489422 | exact (earlier) | 12.7397 | 11.7786 | 1.2503 | 0.2892 |
| 4 | (1,1) | 17.094282 | exact (earlier) | 18.4945 | 17.5236 | 1.4002 | 0.4294 |
| 5 | (3,0) | 21.376199 | exact (earlier) | 22.5996 | 21.6381 | 1.2234 | 0.2619 |
| 6 | (2,1) | 27.697247 | exact (this map) | 29.0197 | 28.0669 | 1.3224 | 0.3696 |
| 7 | (0,2) | 29.639390 | exact (earlier) | 31.1678 | 30.1456 | 1.5284 | 0.5062 |
| 8 | (4,0) | 33.914349 | exact (earlier) | 35.1246 | 34.1616 | 1.2103 | 0.2472 |
| 9 | (1,2) | 39.308624 | exact (this map) | 40.7554 | 39.7665 | 1.4468 | 0.4579 |
| 10 | (3,1) | 41.074347 | exact (this map) | 42.3557 | 41.4085 | 1.2813 | 0.3341 |
| 11 | (5,0) | 49.115375 | exact (this map) | 50.3181 | 49.3534 | 1.2027 | 0.2380 |
| 12 | (2,2) | 52.589164 | exact (this map) | 53.9687 | 53.0027 | 1.3795 | 0.4136 |
| 13 | (4,1) | 57.164950 | exact (this map) | 58.4220 | 57.4757 | 1.2570 | 0.3108 |
| 14 | (0,3) | 59.318660 | exact (earlier) | 60.8467 | 59.8157 | 1.5280 | 0.4970 |
| 15 | (6,0) | 66.983637 | exact (this map) | 68.1816 | 67.2154 | 1.1979 | 0.2318 |
| 16 | (3,2) | 69.013889 | exact (this map) | 70.3474 | 69.3939 | 1.3335 | 0.3800 |
| 17 | (1,3) | 71.191232 | exact (this map) | 72.6642 | 71.6612 | 1.4730 | 0.4700 |
| 18 | (5,1) | 75.948653 | exact (this map) | 77.1899 | 76.2429 | 1.2413 | 0.2942 |
| 19 | (2,3) | 86.787330 | exact (this map) | 88.2061 | 87.2254 | 1.4188 | 0.4381 |
| 20 | (7,0) | 87.521469 | exact (this map) | 88.7160 | 87.7486 | 1.1945 | 0.2271 |
| 21 | (4,2) | 88.354194 | exact (this map) | 89.6564 | 88.7090 | 1.3022 | 0.3548 |
| 22 | (6,1) | 97.416801 | exact (this map) | 98.6472 | 97.6986 | 1.2304 | 0.2818 |
| 23 | (0,4) | 98.265500 | exact (earlier) | 99.7974 | 98.7585 | 1.5319 | 0.4930 |

Order of appearance:

* exact:  (1,0) < (0,1) < (2,0) < (1,1) < (3,0) < (2,1) < (0,2) < (4,0) < (1,2) < (3,1) < (5,0) < (2,2) < (4,1) < (0,3) < (6,0) < (3,2) < (1,3) < (5,1) < (2,3) < (7,0) < (4,2) < (6,1) < (0,4)
* BO+ad: (1,0) < (0,1) < (2,0) < (1,1) < (3,0) < (2,1) < (0,2) < (4,0) < (1,2) < (3,1) < (5,0) < (2,2) < (4,1) < (0,3) < (6,0) < (3,2) < (1,3) < (5,1) < (2,3) < (7,0) < (4,2) < (6,1) < (0,4)
* BO:    (1,0) < (0,1) < (2,0) < (1,1) < (3,0) < (2,1) < (0,2) < (4,0) < (1,2) < (3,1) < (5,0) < (2,2) < (4,1) < (0,3) < (6,0) < (3,2) < (1,3) < (5,1) < (2,3) < (7,0) < (4,2) < (6,1) < (0,4)
* Frolov 1992 II (S, P, D only): (1,0) < (0,1) < (2,0) < (1,1) < (2,1) < (0,2) < (1,2) < (2,2) < (0,3) < (1,3)
* exact restricted to Frolov's levels: (1,0) < (0,1) < (2,0) < (1,1) < (2,1) < (0,2) < (1,2) < (2,2) < (0,3) < (1,3)

## Number of bound levels n(M) of the 1s sigma_g manifold

Exact counting uses the exact M_c (earlier + this map); levels without an exact value yet use the frozen prediction and are counted in the 'of which predicted' column. Valid for M <= 100 (levels with BO+ad M_c in (101, 105] are outside the map).

| M | n exact | of which predicted | n BO | n BO+ad |
|---|---|---|---|---|
| 2 | 1 | 0 | 1 | 1 |
| 5 | 2 | 0 | 1 | 2 |
| 8.88 | 2 | 0 | 2 | 2 |
| 10 | 3 | 0 | 2 | 2 |
| 15 | 4 | 0 | 4 | 4 |
| 17.752 | 5 | 0 | 4 | 5 |
| 20 | 5 | 0 | 5 | 5 |
| 25 | 6 | 0 | 6 | 6 |
| 26.585 | 6 | 0 | 6 | 6 |
| 30 | 8 | 0 | 7 | 7 |
| 40 | 10 | 0 | 9 | 10 |
| 50 | 12 | 0 | 11 | 12 |
| 60 | 15 | 0 | 14 | 15 |
| 70 | 17 | 0 | 16 | 17 |
| 80 | 19 | 0 | 19 | 19 |
| 90 | 22 | 0 | 22 | 22 |
| 100 | 24 | 0 | 24 | 24 |

## Muonic consistency check

Literature bound-level sets (verified against Frolov & Wardlaw 2011 and Frolov 2012): pp mu 2 levels (0,0),(1,0); dd mu 5 levels (+ (0,1),(1,1),(2,0)); tt mu 6 levels (+ (3,0)).

| ion | M | bound levels from the map | literature | agree |
|---|---|---|---|---|
| pp mu | 8.88 | (0,0), (1,0) | (0,0), (1,0) | yes (near: (0,1) M_c=9.526) |
| dd mu | 17.752 | (0,0), (0,1), (1,0), (1,1), (2,0) | (0,0), (0,1), (1,0), (1,1), (2,0) | yes (near: (1,1) M_c=17.094) |
| tt mu | 26.585 | (0,0), (0,1), (1,0), (1,1), (2,0), (3,0) | (0,0), (0,1), (1,0), (1,1), (2,0), (3,0) | yes |

## Orthogonalization check of run 5 (2026-11-17)

Single-pass canonical orthogonalization caused spurious levels for L = 7 at tol <= 1e-14 (threshold_map_L7hp). Run 5 uses tol 1e-12; for three representative standard largest bases the crossing was recomputed with a second orthogonalization pass (check_run5_orthog.py). The single-pass max|X^T S X - I| (2e-3..6e-3) is dominated by rounding in forming X^T S X in double and is not itself an error of the result; the crossing shifts are:

* J2v1 'compact N=16 + 3 diffuse' size=6216 kept=4061: single-pass max|X^T S X - I| = 2.4e-03 (eig 0.997533..1.002484); M* single 27.6972470555 (production 27.6972470555), two-pass 27.6972470555, shift +1.14e-11
* J5v0 'compact N=12 + 2 diffuse' size=5430 kept=3855: single-pass max|X^T S X - I| = 5.5e-03 (eig 0.991572..1.006026); M* single 49.1153803378 (production 49.1153803378), two-pass 49.1153803378, shift +7.60e-11
* J6v1 'compact N=10 + 2 diffuse' size=5638 kept=3531: single-pass max|X^T S X - I| = 4.3e-03 (eig 0.993496..1.005526); M* single 97.4167766625 (production 97.4167766613), two-pass 97.4167766620, shift -4.28e-10

All shifts (<= 4.3e-10) are far below the quoted errors of these targets (1e-9, 4.9e-6, 2.4e-5): run 5 is not affected.
