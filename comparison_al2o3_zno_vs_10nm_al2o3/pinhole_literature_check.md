# Pinhole and literature check for 10 nm Al2O3 reference

## Local model diagnosis

The previously generated 10 nm Al2O3 reference is an extreme "through-pinhole"
case, not a general representation of dense amorphous Al2O3.

Connectivity check of `/home/xiao/文档/KMC/kmc_al2o3_10nm_200k/al2o3_10nm_raw.npz`:

| item | value |
|---|---:|
| grid size | 20 x 640 cells |
| dx | 0.5 nm |
| film thickness | 10 nm |
| permeable pinhole cell fraction | 0.02039 |
| top pinhole mask fraction | 0.028125 |
| connected pinhole components | 4 |
| top-to-bottom percolating pinhole components | 4 |
| Ptrans | 1.0 |
| mean FPT | 0.656 s |

Therefore the earlier bilayer advantage of about 6.2e3 lower WVTR compares the
bilayer against a 10 nm Al2O3 film with four fully connected top-to-bottom
pinhole paths. It demonstrates the danger of direct pinhole leakage, but it
overstates the advantage over a high-quality dense 10 nm Al2O3 film.

Connectivity check of the bilayer result:

| item | value |
|---|---:|
| grid size | 21 x 640 cells |
| permeable cell fraction | 0.19635 |
| connected permeable components | 31 |
| top-to-bottom percolating permeable components | 1 |
| Ptrans | 0.1833 |
| mean FPT | 4258.36 s |

The bilayer has a geometrically connected permeable network, but it is
lateral, tortuous, partly blocked, and has higher ZnO grain-boundary/interface
barriers. The model advantage comes from converting direct Al2O3 pinhole flow
into ZnO grain-boundary-network search and delayed first passage.

## Literature points

1. The 2022 Thin Solid Films Al2O3/ZnO paper assumes dense Al2O3 and ZnO are
   impermeable because the water molecule diameter is larger than the average
   bond lengths; water transport is assigned to defects in amorphous layers and
   grain boundaries in polycrystalline layers.
   Source: https://www.sciencedirect.com/science/article/pii/S0040609022004862

2. Thickness-dependent ALD alumina studies report a critical continuous-film
   thickness around 5-10 nm, and above that range the barrier performance is
   dominated by defect density. For 15-25 nm, both defect density and WVTR can
   decrease exponentially with thickness.

3. UV-ALD Al2O3 studies show that process quality matters strongly. A 10 nm
   UV-ALD Al2O3 film on PET had WVTR around 1.42e-3 g m^-2 day^-1, while a
   10 nm thermal ALD film had WVTR around 6.11e-1 g m^-2 day^-1 under MOCON
   testing. The same paper attributes rapid Ca oxidation in poorer films to
   pinholes forming direct water-vapor pathways.
   Source: https://pubs.rsc.org/en/content/articlehtml/2017/ra/c6ra27759d

4. Plasma-assisted ALD work reported WVTR about 5e-3 g m^-2 day^-1 for a
   20 nm Al2O3 film on PEN, showing that realistic Al2O3 reference values are
   process-dependent and can be orders of magnitude below a deliberately
   through-pinhole model.
   Source: https://research.tue.nl/en/publications/plasma-assisted-atomic-layer-deposition-of-al2o3-moisture-permeat/

## Recommendation for the next comparison

Do not use only the through-pinhole 10 nm Al2O3 model as the main reference.
Use at least three Al2O3 references:

1. Dense/high-quality 10 nm Al2O3 calibrated to literature WVTR.
2. Defective but non-through or partially staggered 10 nm Al2O3.
3. Worst-case through-pinhole 10 nm Al2O3, clearly labeled as an upper-leakage
   bound.

Then compare the bilayer against this uncertainty envelope rather than a
single extreme baseline.
