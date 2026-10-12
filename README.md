# Multi-Scale Intersection for Improved Edge Separation in Robust SAR Edge Detection

Aditya Rane , Jay Gandhi 
Mukesh Patel School of Technology Management and Engineering, NMIMS, Mumbai, India

---

## Overview

ROSED detects edges in SAR images with statistical tests, but it uses one fixed
14 x 14 window and, like earlier ratio detectors, assumes a single edge inside
that window. Two boundaries closer together than the window are merged into one
response.

This work reimplements ROSED, verifies each component against theory and
reproduces its published results. It then computes the Ratio of Medians at three
window sizes, s = 8, 14 and 20, and keeps only the detections that all three
agree on, before the likelihood ratio test localises the edge. The test emits one
edge per detection run, so requiring agreement holds the runs to the width of the
smallest window and two close boundaries stay in separate runs.

The only change from ROSED is where the scales are combined: the logical AND is
applied to the ROM maps before run grouping and the LRT, not after.

## Result

On simulated images with two boundaries 10 pixels apart, the intersection
resolves both boundaries in 100.0 +/- 0.1% of image rows, against 24.1 +/- 3.3%
for the original method, over ten speckle realisations. It gives 18.7 to 23.6%
fewer false positives than the best single window at every separation tested, in
70 of 70 cases.

Two further findings are reported: two evaluation parameters in the original
paper contradict each other, and on open 8-bit SAR products 25 to 37% of pixels
sit at the minimum value, which disables the finest window and makes the CFAR
false alarm rate vary about fifty-fold with brightness.

Full results are in the paper.

---

## Repository layout

```
code/      everything that runs
data/      the dataset
results/   the frozen tables, figures and logs

README.md
```

### code/

The shared library is a single file. Every function that more than one
experiment uses lives in `rosed.m` as a static method, so all results come from
one copy of the code. A caller reaches it with one line, `import rosed.*`, after
which every function is called by its plain name. Each local function inside a
script needs its own copy of that line, because the import is scoped to one
function.

| | |
| --- | --- |
| `rosed.m` | the shared library: 37 static methods in seven chapters, covering speckle simulation and the test images, the detector, the scale combination, calibration, the evaluation metrics, the regression gate, and small utilities |
| `run_all_simulated.m` | runs the seven simulated experiments in order and logs each one |
| `rosed_pipeline.m` | the reimplementation end to end, with its own regression check |
| `fig6e_extended.m` | image E across separations of 6 to 24 pixels, plus run widths |
| `fig6e_robustness.m` | image E over ten speckle realisations |
| `scale_sets_experiment.m` | image E with other window-size sets |
| `fig6c_guardcells.m` | the dark-line test with guard cells excluded |
| `classical_baselines.m` | Sobel, Laplacian of Gaussian and Canny |
| `morphology_variant.m` | closing and small-segment removal |
| `run_pooled_calibration.m` | real-data thresholds pooled over the calibration regions |
| `real_calibrate_neg.m` | the negative CFAR constant, which the pooled run does not produce |
| `real_detect.m` | every strategy on the three real crops |
| `real_channels.m` | the ROM and CFAR channels separately, and the CFAR rate by brightness |
| `diagnose_patches.m` | per-region thresholds, used to decide whether pooling was valid |
| `test_ordfilt2.m` | proves the toolbox order-statistic filter matches the hand-written one |
| `figures_code/` | the seven scripts that draw the figures in the paper |
| `extract_crops.py` | cuts the crops listed in `data/crops.csv` and writes the previews |
| `find_clean_patches.py` | searches the scenes for homogeneous calibration regions |
| `auto_pick_crops.py` | ranks candidate test crops by vertical structure |
| `refine_calib.py` | refines a calibration region |
| `results_table.csv`, `results_fig6d.csv`, `results_fig6e.csv` | the regression gate's baseline. Do not remove them. A missing baseline is treated as a failure, not as a skip |

### results/

| | |
| --- | --- |
| `tables/` | 19 CSV files and 4 `.mat` threshold files. Every number reported in the paper comes from one of these |
| `figures/` | `paper/` holds the figures used in the paper, `real/` the real-data output, and the loose PNGs the per-experiment plots |
| `logs/` | the console output of each run, one file per experiment |

---

## The dataset

```
data/
  crops.csv     the crop definitions: which scene, which region, what label
  crops/        12 calibration regions, 3 test crops, and excluded/ with the 3
                regions trimmed from the pool
  previews/     one downsampled preview per scene, for picking crops by eye
  simulated/    the saved simulated test images, so the figures can be redrawn
                without generating new random numbers
  raw/          the nine source scenes. NOT included, see below
```

The real images come from **Umbra Open Data**, released under **CC BY 4.0**:

> https://umbra.space/open-data

Nine scenes in the geocoded GEC GeoTIFF format are used, from 8001 x 8001 to
22859 x 22859 pixels, covering cities, a port, a harbour and farmland.

### Getting the scenes

`data/raw/` is not in this repository. The nine scenes come to 3.08 GB and six
of them are individually larger than the 100 MB limit for a single file. They are
public, so downloading them is better than re-hosting them.

Download these and put them in `data/raw/`, keeping the names exactly:

```
2025-05-17-17-43-46_UMBRA-08_GEC.tif
2025-06-15-08-28-22_UMBRA-08_GEC.tif
2025-06-22-23-57-52_UMBRA-10_GEC.tif
2025-07-25-05-15-29_UMBRA-07_GEC.tif
2025-10-13-05-31-39_UMBRA-07_GEC.tif
2025-10-31-06-42-11_UMBRA-10_GEC.tif
2025-12-18-18-07-38_UMBRA-07_GEC.tif
2025-12-19-12-49-51_UMBRA-07_GEC.tif
2026-01-02-03-22-31_UMBRA-09_GEC.tif
```

Everything derived from them is included, so the work reproduces without the
3 GB. `data/previews/` shows what each scene looks like and `data/crops/` holds
the regions every real-data result was computed from. The scenes are needed only
to cut fresh crops.

---

## Pipeline

```
SIMULATED BRANCH

   rosed.m                       speckle from the physical model, fixed seeds
      |
      v
   data/simulated/*.mat          saved test images A to E
      |                          run_all_simulated
      v
   results/tables/*.csv          frozen result tables
   results/figures/*.png
   results/logs/*.txt


REAL BRANCH

   data/raw/*.tif                9 Umbra scenes  (download separately)
      |                          extract_crops.py, driven by data/crops.csv
      v
   data/crops/*.tif              12 calibration regions, 3 test crops
      |                          run_pooled_calibration, real_calibrate_neg
      v
   results/tables/calibration_real_pooled.mat     thresholds
      |                          real_detect, real_channels
      v
   results/tables/real_*.csv
   results/figures/real/


FIGURES FOR THE PAPER

   data/simulated/*.mat
   data/crops/*real_pair*.tif
      |                          figures_code/fig_test_images.m
      |                          figures_code/paper_figures.m
      v
   results/figures/paper/*.png
```

### Running it

**Simulated.** From MATLAB, in the `code` folder:

```matlab
cd code
run_all_simulated
```

About one and a half minutes. Three regression checks must print:

```
[regression] PASS: 12 match, 0 differ, 0 missing   results_table.csv
[regression] PASS:  4 match, 0 differ, 0 missing   results_fig6d.csv
[regression] PASS: 16 match, 0 differ, 0 missing   results_fig6e.csv
```

If any fails, stop. It means a published number no longer reproduces.

**Real data.** The crops are already cut, so the first two steps are needed only
if the scenes are re-downloaded:

```
python code/extract_crops.py survey     writes the previews
python code/extract_crops.py extract    cuts the crops in data/crops.csv
```

then in MATLAB from `code`, in this order:

```matlab
run_pooled_calibration     ROM thresholds and the positive CFAR constants
real_calibrate_neg         the negative CFAR constant
real_detect                every strategy on the three test crops
real_channels              the ROM and CFAR channels separately
```

**Figures.** From MATLAB in `code`:

```matlab
addpath code/figures_code
fig_test_images
paper_figures
```

These print six numbers as they run. They must read 23.4% with 228 false
positives, 99.6% with 164, then 3.81% and 0.69%, then 2.66, 20.18 and 30.61%. If
any differs, the figures and the tables are coming from different images.

---

## Verification

The results are frozen. Three things enforce it.

Each experiment compares its own output against a saved baseline every time it
runs, and stops if anything differs, so the code cannot be reorganised and
quietly change a published number. A missing baseline counts as a failure rather
than a skip, so a check that cannot run cannot pass.

Every number in the paper is checked against the CSV files in `results/tables/`
at the paper's own precision, 292 checks in total.

Every random seed is a function of the loop indices alone, so the repeated
realisations give identical results whether or not the Parallel Computing
Toolbox is present.

The real scenes have no ground truth, since nobody has marked where the true
edges are, so those results are reported qualitatively together with
measurements that do not need ground truth.

---

## Requirements

MATLAB R2025a or later with the Image Processing, Statistics and Machine
Learning, and Parallel Computing toolboxes. The parallel toolbox is optional.

Python 3.13 with `numpy`, `tifffile`, `imagecodecs` and `pillow`, used only to
read the GeoTIFF scenes and cut the crops.

---

## Reference

The method extended here is:

> A. Hashem et al., "Robust SAR Edge Detection", European Radar Conference
> (EuRAD), 2025.

The full reference list is in the paper.

## Licence

The Umbra scenes are CC BY 4.0, Umbra Lab Inc. The code in this repository is
released for academic use.
