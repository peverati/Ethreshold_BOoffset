# Result: exact critical masses of L > 0 levels of {X+, X+, e-} (thresholds_Lv)

*Last updated: 2026-11-15 22:45*

Written by `analyze.py` from the run_exact.py CSV files. Crossings M* of fixed bases are rigorous upper bounds to M_c (variational root below E_th proves binding), up to rounding.

## Blind-test integrity

* sha256(PREDICTION_J1v1.md) = 551b8f2899a56b9056eaac6a42ea8d6b48e78eb22fcfc6235d156ddb5853e887: matches prediction_J1v1.sha256.
* Frozen 2026-11-15 20:17. Caveat (PROGRESS.md): the J1v1 scale optimization (stage 1, N = 12) ran 20:11-20:13, before the freeze; it computed no crossing.

## (J=1, v=1): L=1 'A', root 1

| basis | size | kept | M* (upper bound) |
|---|---|---|---|
| compact N=12 | 910 | 860 | 17.1488526287 |
| compact N=14 | 1360 | 1264 | 17.1287130315 |
| compact N=16 | 1938 | 1777 | 17.1173197674 |
| compact N=18 | 2660 | 2405 | 17.1104266503 |
| compact N=20 | 3542 | 3168 | 17.1061007438 |
| compact N=18 + 1 diffuse | 3050 | 2570 | 17.0984432828 |
| compact N=18 + 2 diffuse | 3560 | 2845 | 17.0944235358 |
| compact N=18 + 3 diffuse | 4274 | 3313 | 17.0942872329 |
| compact N=18 + 4 diffuse | 5156 | 3756 | 17.0942820979 |
| compact N=18 + 5 diffuse | 6332 | 4445 | 17.0942819121 |
| compact N=16 + 5 diffuse | 5610 | 3856 | 17.0942819121 |
| compact N=20 + 5 diffuse | 7214 | 5166 | 17.0942819120 |

* **best proven upper bound (up to rounding): 17.094281912**
* rounding-safe bound (margin 1e-11 Eh): 17.094281939
* tol 1e-10 on 'compact N=18 + 5 diffuse': 17.0942819137 (difference from tol 1e-12 bound: +1.8e-09)
* compact-size sequence with 5 diffuse sets, N = [16, 18, 20]: 17.094281912, 17.094281912, 17.094281912; ratio nan
* **extrapolated M_c = 17.09428191 +/- 1.0e-09** (geometric in N; assumes geometric convergence, the bound above is rigorous)
* diffuse increments at N_D: 1:17.0984432828, 2:17.0944235358, 3:17.0942872329, 4:17.0942820979, 5:17.0942819121; last increment -1.9e-07

Binding energies (largest basis):

| M | x = M/M* - 1 | E_th - E | E_b/(M - M*) |
|---|---|---|---|
| 17.10282905 | 0.0005 | 3.322414e-06 | 3.8872e-04 |
| 17.11137619 | 0.0010 | 6.769209e-06 | 3.9599e-04 |
| 17.12847048 | 0.0020 | 1.389492e-05 | 4.0642e-04 |
| 17.17975332 | 0.0050 | 3.654181e-05 | 4.2753e-04 |
| 17.26522473 | 0.0100 | 7.724617e-05 | 4.5188e-04 |
| 17.43616755 | 0.0200 | 1.665072e-04 | 4.8703e-04 |
| 17.94899601 | 0.0500 | 4.764764e-04 | 5.5747e-04 |
| 18.80371010 | 0.1000 | 1.083833e-03 | 6.3403e-04 |
| 20.51313829 | 0.2000 | 2.489621e-03 | 7.2820e-04 |

* quadratic fit (x <= 0.05) root 17.09700359; slope at root 4.4124e-04 Eh per mass unit (linear threshold law for L > 0)

* BO 18.494475, BO+ad 17.523638
* literature / consistency: Frolov & Wardlaw arXiv:1008.3010: bound up to m_mu ~ 214.72 m_e, M_c ~ 17.09 (basis-limited upper estimate)
* non-adiabatic offset M_c(BO+ad) - M_c(exact) = 0.42936 (predicted 0.41647); offset vs BO 1.40019

## (J=3, v=0): L=3 'A', root 0

| basis | size | kept | M* (upper bound) |
|---|---|---|---|
| compact N=8 | 660 | 608 | 21.3764296019 |
| compact N=10 | 1144 | 1013 | 21.3762394874 |
| compact N=12 | 1820 | 1565 | 21.3762077927 |
| compact N=14 | 2720 | 2278 | 21.3762014224 |
| compact N=12 + 1 diffuse | 2600 | 2012 | 21.3762033839 |
| compact N=12 + 2 diffuse | 3620 | 2578 | 21.3762034641 |
| compact N=12 + 3 diffuse | 5048 | 3496 | 21.3762033518 |
| compact N=12 + 4 diffuse | 6812 | 4356 | 21.3762033793 |
| compact N=10 + 4 diffuse | 6136 | 3897 | 21.3762219825 |
| compact N=14 + 4 diffuse | 7712 | 4982 | 21.3761998653 |

* **best proven upper bound (up to rounding): 21.376199865**
* rounding-safe bound (margin 1e-11 Eh): 21.376203385
* tol 1e-10 on 'compact N=12 + 4 diffuse': 21.3762335840 (difference from tol 1e-12 bound: +3.4e-05)
* compact-size sequence with 4 diffuse sets, N = [10, 12, 14]: 21.376221983, 21.376203379, 21.376199865; ratio 0.189
* **extrapolated M_c = 21.37619905 +/- 8.2e-07** (geometric in N; assumes geometric convergence, the bound above is rigorous)
* diffuse increments at N_D: 1:21.3762033839, 2:21.3762034641, 3:21.3762033518, 4:21.3762033793; last increment +2.7e-08

Binding energies (largest basis):

| M | x = M/M* - 1 | E_th - E | E_b/(M - M*) |
|---|---|---|---|
| 21.38688797 | 0.0005 | 1.829536e-05 | 1.7118e-03 |
| 21.39757607 | 0.0010 | 3.660069e-05 | 1.7122e-03 |
| 21.41895226 | 0.0020 | 7.322263e-05 | 1.7127e-03 |
| 21.48308086 | 0.0050 | 1.831700e-04 | 1.7138e-03 |
| 21.58996186 | 0.0100 | 3.666386e-04 | 1.7152e-03 |
| 21.80372386 | 0.0200 | 7.341244e-04 | 1.7172e-03 |
| 22.44500986 | 0.0500 | 1.837042e-03 | 1.7188e-03 |
| 23.51381985 | 0.1000 | 3.660902e-03 | 1.7126e-03 |
| 25.65143984 | 0.2000 | 7.195870e-03 | 1.6832e-03 |

* quadratic fit (x <= 0.05) root 21.37628244; slope at root 1.7160e-03 Eh per mass unit (linear threshold law for L > 0)

* BO 22.599590, BO+ad 21.638056
* literature / consistency: ddmu (M = 17.752) has no bound J=3 level, ttmu (26.585) has one: 17.752 < M_c < 26.585
* non-adiabatic offset M_c(BO+ad) - M_c(exact) = 0.26186 (predicted 0.25622); offset vs BO 1.22339

## (J=4, v=0): L=4 'S', root 0

| basis | size | kept | M* (upper bound) |
|---|---|---|---|
| compact N=8 | 850 | 767 | 33.9144342591 |
| compact N=10 | 1466 | 1270 | 33.9143594898 |
| compact N=12 | 2324 | 1937 | 33.9143502790 |
| compact N=10 + 1 diffuse | 2636 | 1917 | 33.9143589247 |
| compact N=10 + 2 diffuse | 4166 | 2625 | 33.9143588821 |
| compact N=10 + 3 diffuse | 6308 | 3767 | 33.9143596079 |
| compact N=8 + 3 diffuse | 5692 | 3362 | 33.9144200567 |
| compact N=12 + 3 diffuse | 7166 | 4339 | 33.9143502902 |

* **best proven upper bound (up to rounding): 33.914350279**
* rounding-safe bound (margin 1e-11 Eh): 33.914359616
* tol 1e-10 on 'compact N=10 + 3 diffuse': 33.9144855202 (difference from tol 1e-12 bound: +1.4e-04)
* compact-size sequence with 3 diffuse sets, N = [8, 10, 12]: 33.914420057, 33.914359608, 33.914350290; ratio 0.154
* **extrapolated M_c = 33.91434859 +/- 1.7e-06** (geometric in N; assumes geometric convergence, the bound above is rigorous)
* diffuse increments at N_D: 1:33.9143589247, 2:33.9143588821, 3:33.9143596079; last increment +7.3e-07

Binding energies (largest basis):

| M | x = M/M* - 1 | E_th - E | E_b/(M - M*) |
|---|---|---|---|
| 33.93130745 | 0.0005 | 2.053588e-05 | 1.2110e-03 |
| 33.94826463 | 0.0010 | 4.108204e-05 | 1.2113e-03 |
| 33.98217898 | 0.0020 | 8.217109e-05 | 1.2115e-03 |
| 34.08392203 | 0.0050 | 2.054082e-04 | 1.2113e-03 |
| 34.25349378 | 0.0100 | 4.106803e-04 | 1.2109e-03 |
| 34.59263728 | 0.0200 | 8.206099e-04 | 1.2098e-03 |
| 35.61006779 | 0.0500 | 2.043140e-03 | 1.2049e-03 |
| 37.30578531 | 0.1000 | 4.046162e-03 | 1.1931e-03 |
| 40.69722033 | 0.2000 | 7.883568e-03 | 1.1623e-03 |

* quadratic fit (x <= 0.05) root 33.91442352; slope at root 1.2131e-03 Eh per mass unit (linear threshold law for L > 0)

* BO 35.124626, BO+ad 34.161564
* non-adiabatic offset M_c(BO+ad) - M_c(exact) = 0.24722 (predicted 0.24406); offset vs BO 1.21028

## 2p sigma_u J=1: L=1 'S', root 0

| basis | size | kept | M* (upper bound) |
|---|---|---|---|
| compact N=8 | 330 | 316 | nan |
| compact N=10 | 572 | 531 | nan |
| compact N=10 + 1 diffuse | 1802 | 1496 | 815.9867892238 |
| compact N=10 + 2 diffuse | 3032 | 1910 | 813.9078470426 |
| compact N=10 + 3 diffuse | 4562 | 2471 | 813.6540524612 |
| compact N=10 + 4 diffuse | 6092 | 2909 | 813.6488581679 |
| compact N=10 + 5 diffuse | 7922 | 3513 | 813.6470878126 |
| compact N=8 + 5 diffuse | 7680 | 3362 | 813.7250061031 |
| compact N=12 + 5 diffuse | 8260 | 3737 | 813.6015319265 |

* **best proven upper bound (up to rounding): 813.601531926**
* rounding-safe bound (margin 1e-11 Eh): 813.648514641
* tol 1e-10 on 'compact N=10 + 5 diffuse': 813.7006388165 (difference from tol 1e-12 bound: +9.9e-02)
* compact-size sequence with 5 diffuse sets, N = [8, 10, 12]: 813.725006103, 813.647087813, 813.601531926; ratio 0.585
* **extrapolated M_c = 813.53740385 +/- 6.4e-02** (geometric in N; assumes geometric convergence, the bound above is rigorous)
* diffuse increments at N_D: 1:815.9867892238, 2:813.9078470426, 3:813.6540524612, 4:813.6488581679, 5:813.6470878126; last increment -1.8e-03

Binding energies (largest basis):

| M | x = M/M* - 1 | E_th - E | E_b/(M - M*) |
|---|---|---|---|
| 814.00833269 | 0.0005 | 2.591989e-09 | 6.3716e-09 |
| 814.41513346 | 0.0010 | 5.591212e-09 | 6.8722e-09 |
| 815.22873499 | 0.0020 | 1.174696e-08 | 7.2191e-09 |
| 817.66953959 | 0.0050 | 3.104824e-08 | 7.6323e-09 |
| 821.73754725 | 0.0100 | 6.510217e-08 | 8.0017e-09 |
| 829.87356256 | 0.0200 | 1.380228e-07 | 8.4822e-09 |
| 854.28160852 | 0.0500 | 3.810952e-07 | 9.3681e-09 |
| 894.96168512 | 0.1000 | 8.338956e-07 | 1.0249e-08 |
| 976.32183831 | 0.2000 | 1.824951e-06 | 1.1215e-08 |

* quadratic fit (x <= 0.05) root 813.74736165; slope at root 7.9366e-09 Eh per mass unit (linear threshold law for L > 0)

* BO 815.558715 (BO+ad not available for 2p sigma_u)
* literature / consistency: none checked

## Summary

| target | BO | BO+ad | best upper bound | extrapolated M_c | prediction (frozen) | exact - pred | pass |
|---|---|---|---|---|---|---|---|
| J1v1 | 18.49448 | 17.52364 | 17.0942819 | 17.0942819 | 17.1072 +/- 0.0295 | -0.0129 | PASS |
| L3A | 22.59959 | 21.63806 | 21.3761999 | 21.3761990 | 21.3818 +/- 0.0295 | -0.0056 | PASS |
| L4S | 35.12463 | 34.16156 | 33.9143503 | 33.9143486 | 33.9175 +/- 0.0295 | -0.0032 | PASS |
| L1S | 815.55871 | None | 813.6015319 | 813.5374038 | | | |
