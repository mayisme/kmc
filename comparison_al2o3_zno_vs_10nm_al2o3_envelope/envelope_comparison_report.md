# Bilayer vs Al2O3 reference envelope

This comparison replaces the earlier single extreme through-pinhole baseline with an envelope of Al2O3 references.

| Case | Type | Ptrans | FPT mean (s) | tau mean | D_eff (m2/s) | WVTR C1 (g m-2 day-1) |
|---|---|---:|---:|---:|---:|---:|
| 4.5 nm Al2O3 / 6 nm ZnO bilayer KMC | KMC result | 0.1833 | 4258.36 | 1.147e+05 | 1.295e-20 | 7.211e-04 |
| 10 nm high-quality UV-ALD Al2O3 literature reference | literature-calibrated | N/A | N/A | N/A | 2.428e-20 | 0.00142 |
| 20 nm PA-ALD Al2O3 literature context | literature-calibrated context | N/A | N/A | N/A | 8.548e-20 | 0.005 |
| 10 nm thermal ALD Al2O3 poorer-process reference | literature-calibrated | N/A | N/A | N/A | 1.045e-17 | 0.611 |
| 10 nm Al2O3 staggered non-through defect idealization | connectivity-limited KMC geometry | 0 | N/A | N/A | 0 | 0 |
| 10 nm Al2O3 through-pinhole worst-case KMC | KMC extreme leakage bound | 1 | 0.655689 | 331.414 | 7.626e-17 | 4.46041 |

## Relations to bilayer

| Reference | WVTR relation | D_eff relation |
|---|---:|---:|
| 10 nm high-quality UV-ALD Al2O3 literature reference | 1.97x lower | 1.88x lower |
| 20 nm PA-ALD Al2O3 literature context | 6.93x lower | 6.6x lower |
| 10 nm thermal ALD Al2O3 poorer-process reference | 847x lower | 807x lower |
| 10 nm Al2O3 staggered non-through defect idealization | bilayer higher than ideal zero-leakage reference | bilayer higher than ideal zero-leakage reference |
| 10 nm Al2O3 through-pinhole worst-case KMC | 6.19e+03x lower | 5.89e+03x lower |

## Interpretation

- The previous 6.2e3 advantage is valid only against the through-pinhole worst-case Al2O3 model.
- Against a high-quality 10 nm UV-ALD literature reference, the current bilayer KMC WVTR is about 2x lower, not thousands of times lower.
- Against a poorer 10 nm thermal ALD literature reference, the bilayer is hundreds of times lower.
- Against a perfect non-through-defect idealization, any finite bilayer leakage is higher than the zero-leakage reference.
- Therefore the defensible claim is an envelope: bilayer performance is strongly better than through-defect or poorer-process Al2O3, comparable to or moderately better than high-quality 10 nm UV-ALD Al2O3 under the chosen calibration, but not better than an ideal fully non-percolating dense film.
