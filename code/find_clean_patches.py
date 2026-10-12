#!/usr/bin/env python3
"""
find_clean_patches.py — locate uniform LAND patches suitable for calibration

Run from the project root:
    python code/find_clean_patches.py

FILTER HISTORY
--------------
v1  scored on closeness to CV = 0.523. Selected the maritime scene
    exclusively, because ocean is the only genuinely Rayleigh surface
    in this collection. Calibrating on water and detecting on land is
    not defensible.

v2  excluded the maritime scene and scored on quadrant uniformity
    instead. Produced a good spread across 8 land scenes, but the
    per-patch diagnostic showed tau varying 43x across the pool, so
    pooling was invalid. Cause: every patch has 18-39% of its pixels at
    value <= 3. In uint8 data with a mean near 30, a half-window median
    can land on 2 while the other half sits at 25, giving r-hat = 0.08
    from quantization alone with no edge present. Patches above roughly
    23% low-value pixels had their tails set by this effect, collapsing
    tau by an order of magnitude.

v3  (this version) adds LOW_FRAC_MAX. Patches whose low-value fraction
    exceeds the limit are rejected, because their r-hat tail measures
    quantization rather than clutter. The five patches that passed this
    criterion in the v2 pool had tau spanning only 1.9x - comfortably
    poolable - and sat in the same order of magnitude as the simulated
    reference.

Fewer, consistent patches are worth more than many inconsistent ones.
If the yield is too low, raise LOW_FRAC_MAX in small steps and re-check
the tau ratio with diagnose_patches.m after every change.
"""

import os
import glob
import numpy as np
import tifffile

RAW      = "data/raw"
CROPFILE = "data/crops.csv"

PATCH  = 128
STEP   = 256              # finer stride than v2, to offset the stricter filter
N_KEEP = 24
PER_SCENE_MAX = 6         # raised, since fewer patches will qualify

CV_MIN, CV_MAX  = 0.55, 1.05    # plausible land clutter
SPREAD_MAX      = 1.20          # quadrant CV agreement
MEAN_SPREAD_MAX = 1.15          # quadrant brightness agreement
LOW_FRAC_MAX    = 22.0          # max % of pixels <= 3   <-- the v3 addition
LOW_VALUE       = 3

# maritime scene - excluded, see v1 note
EXCLUDE = ["2025-06-22-23-57-52_UMBRA-10"]


def quadrant_stats(p):
    """Return (cv_spread, mean_spread) across the four quadrants.

    cv_spread   = max/min of quadrant CVs    -> texture consistency
    mean_spread = max/min of quadrant means  -> brightness consistency

    A patch straddling an edge fails the second even when it passes the
    first, so both are checked.
    """
    h = p.shape[0] // 2
    cvs, mus = [], []
    for q in (p[:h, :h], p[:h, h:], p[h:, :h], p[h:, h:]):
        q = q[q > 0]
        if q.size < 100:
            return 99.0, 99.0
        mu = max(q.mean(), 1e-9)
        cvs.append(q.std() / mu)
        mus.append(mu)
    return max(cvs)/max(min(cvs), 1e-9), max(mus)/max(min(mus), 1e-9)


def main():
    files = sorted(glob.glob(os.path.join(RAW, "*_GEC.tif")))
    files = [f for f in files
             if not any(x in os.path.basename(f) for x in EXCLUDE)]

    if not files:
        print(f"No usable *_GEC.tif in {RAW}/")
        return

    print(f"scanning {len(files)} scenes (maritime excluded)")
    print(f"filters: CV {CV_MIN}-{CV_MAX}, qCV<={SPREAD_MAX}, "
          f"qMean<={MEAN_SPREAD_MAX}, low-frac<={LOW_FRAC_MAX}%\n")

    found = []
    rejected_low = 0

    for f in files:
        name = os.path.basename(f)[:-4]
        print(f"  {name} ...", flush=True)

        with tifffile.TiffFile(f) as t:
            a = t.pages[0].asarray()

        H, W = a.shape
        kept = 0
        rej  = 0
        for r in range(0, H - PATCH, STEP):
            for c in range(0, W - PATCH, STEP):
                p = a[r:r+PATCH, c:c+PATCH].astype(float)
                v = p[p > 0]
                if v.size < 0.98 * PATCH * PATCH:
                    continue

                cv = v.std() / max(v.mean(), 1e-9)
                if cv < CV_MIN or cv > CV_MAX:
                    continue

                # v3: reject patches whose tail would be set by quantization
                lowf = 100.0 * np.count_nonzero(p <= LOW_VALUE) / p.size
                if lowf > LOW_FRAC_MAX:
                    rej += 1
                    continue

                spread, mspread = quadrant_stats(p)
                if spread > SPREAD_MAX or mspread > MEAN_SPREAD_MAX:
                    continue

                score = (spread - 1.0) + 0.5*(mspread - 1.0) + 0.02*lowf
                found.append((score, cv, spread, mspread, lowf,
                              name, c, r, v.mean()))
                kept += 1

        rejected_low += rej
        print(f"      {kept} kept, {rej} rejected on low-value fraction")
        del a

    if not found:
        print(f"\nNothing qualified ({rejected_low} rejected on low-value "
              f"fraction). Raise LOW_FRAC_MAX and retry.")
        return

    found.sort()

    keep, per_scene = [], {}
    for item in found:
        nm = item[5]
        if per_scene.get(nm, 0) >= PER_SCENE_MAX:
            continue
        per_scene[nm] = per_scene.get(nm, 0) + 1
        keep.append(item)
        if len(keep) >= N_KEEP:
            break

    print(f"\n{len(found)} candidates, keeping {len(keep)} "
          f"(max {PER_SCENE_MAX} per scene)\n")
    print(f"{'scene':34s} {'x':>7s} {'y':>7s} {'CV':>6s} {'qCV':>6s} "
          f"{'qMean':>6s} {'low%':>6s} {'mean':>6s}")
    for _, cv, sp, ms, lf, name, c, r, m in keep:
        print(f"{name[:34]:34s} {c:7d} {r:7d} {cv:6.3f} {sp:6.3f} "
              f"{ms:6.3f} {lf:6.2f} {m:6.1f}")

    print("\nscenes represented:")
    for nm, n in sorted(per_scene.items()):
        print(f"  {n} x {nm}")

    lfs = [k[4] for k in keep]
    print(f"\nlow-value fraction in kept set: "
          f"min {min(lfs):.2f}%  max {max(lfs):.2f}%")
    print("If this range is narrow, the pool should be internally consistent.")

    if os.path.exists(CROPFILE):
        text = open(CROPFILE).read()
    else:
        text = "scene,label,x,y,width,height\n"
    lines = [ln for ln in text.splitlines()
             if ln.strip() and ",clean" not in ln]

    for i, (_, cv, sp, ms, lf, name, c, r, m) in enumerate(keep):
        lines.append(f"{name},clean{i:02d},{c},{r},{PATCH},{PATCH}")

    with open(CROPFILE, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    print(f"\n{len(keep)} clean* rows written to {CROPFILE}")
    print("\nNext:")
    print("  del data\\crops\\*clean*.tif")
    print("  python code/extract_crops.py extract")
    print("  then in MATLAB: diagnose_patches('../data/crops','*clean*.tif',14,1e-3)")


if __name__ == "__main__":
    main()
