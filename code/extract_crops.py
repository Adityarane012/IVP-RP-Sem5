#!/usr/bin/env python3
"""
extract_crops.py — Umbra GEC inspection and crop extraction

Two modes:

  python extract_crops.py survey
      Inspects every *_GEC.tif in data/raw/, prints structure and statistics,
      and writes a downsampled preview PNG per scene to data/previews/.
      Use the previews to choose crop coordinates.

  python extract_crops.py extract
      Reads crops.csv and writes each crop to data/crops/ as a plain
      uint8 TIFF that MATLAB's imread opens directly.

crops.csv format (create after surveying):
    scene,label,x,y,width,height
    2025-12-19-12-49-51_UMBRA-07,boundary,3200,4100,512,512
    2025-12-19-12-49-51_UMBRA-07,urban,1800,2400,512,512

Coordinates are in FULL RESOLUTION pixels. The preview PNGs are
downsampled, so multiply preview coordinates by the printed scale factor.
"""

import os
import sys
import csv
import glob

import numpy as np
import tifffile
from PIL import Image

RAW      = "data/raw"
PREVIEWS = "data/previews"
CROPS    = "data/crops"
CROPFILE = "data/crops.csv"

PREVIEW_MAX = 1200          # longest side of preview PNG


# ----------------------------------------------------------------------
def sar_stretch(a):
    """Log scale with 2/98 percentile clipping. DISPLAY ONLY.

    The log destroys the intensity statistics that ROM and OS-CFAR
    assume. Never feed the output of this into a detector.
    """
    a = np.asarray(a, dtype=float)
    L = np.log10(a + 1e-9)
    f = L[np.isfinite(L)]
    lo, hi = np.percentile(f, [2, 98])
    g = np.clip((L - lo) / max(hi - lo, 1e-9), 0, 1)
    return (g * 255).astype(np.uint8)


def homogeneity_scan(a, size=256, step=400):
    """Return the lowest-CV patch found. Useful for picking a region to
    calibrate detection thresholds on real data."""
    best = None
    H, W = a.shape
    for r in range(0, H - size, step):
        for c in range(0, W - size, step):
            p = a[r:r+size, c:c+size].astype(float)
            p = p[p > 0]
            if p.size < 0.9 * size * size:
                continue
            cv = p.std() / max(p.mean(), 1e-9)
            if best is None or cv < best[0]:
                best = (cv, r, c, p.mean())
    return best


# ----------------------------------------------------------------------
def survey():
    os.makedirs(PREVIEWS, exist_ok=True)
    files = sorted(glob.glob(os.path.join(RAW, "*_GEC.tif")))

    if not files:
        print(f"No *_GEC.tif found in {RAW}/")
        return

    print(f"Found {len(files)} scenes\n")

    for f in files:
        name = os.path.basename(f)[:-4]
        size_mb = os.path.getsize(f) / 1e6
        print("=" * 72)
        print(f"{name}   ({size_mb:.0f} MB)")

        with tifffile.TiffFile(f) as t:
            n = len(t.pages)
            p0 = t.pages[0]
            print(f"  full resolution : {p0.shape[0]} x {p0.shape[1]}  dtype={p0.dtype}")
            print(f"  pyramid levels  : {n}")

            # pick a level small enough for a preview
            lvl = 0
            for i, p in enumerate(t.pages):
                if max(p.shape[:2]) <= PREVIEW_MAX * 2:
                    lvl = i
                    break
            else:
                lvl = n - 1
            a = t.pages[lvl].asarray()

        scale = p0.shape[0] / a.shape[0]
        print(f"  preview level   : {lvl}  ({a.shape[0]}x{a.shape[1]}), "
              f"scale factor {scale:.0f}")

        v = a[a > 0].astype(float)
        print(f"  min/mean/max    : {v.min():.0f} / {v.mean():.1f} / {v.max():.0f}")
        print(f"  global CV       : {v.std()/v.mean():.4f}")

        hom = homogeneity_scan(a, size=128, step=100)
        if hom:
            cv, r, c, m = hom
            print(f"  most homogeneous 128x128 (preview coords): "
                  f"CV={cv:.3f} at ({c},{r}) mean={m:.1f}")
            print(f"    -> full-res approx: ({int(c*scale)},{int(r*scale)})")

        out = os.path.join(PREVIEWS, name + ".png")
        Image.fromarray(sar_stretch(a)).save(out)
        print(f"  preview written : {out}")
        print(f"  NOTE: multiply preview coords by {scale:.0f} for crops.csv")

    print("\n" + "=" * 72)
    print("Next: open the previews, pick regions, fill in crops.csv, then run")
    print("  python extract_crops.py extract")


# ----------------------------------------------------------------------
def extract():
    os.makedirs(CROPS, exist_ok=True)

    if not os.path.exists(CROPFILE):
        print(f"{CROPFILE} not found. Create it with header:")
        print("  scene,label,x,y,width,height")
        return

    cache = {}
    written = 0

    with open(CROPFILE) as fh:
        for row in csv.DictReader(fh):
            scene = row["scene"].strip()
            label = row["label"].strip()
            x, y = int(row["x"]), int(row["y"])
            w, h = int(row["width"]), int(row["height"])

            path = os.path.join(RAW, scene + ".tif")
            if not os.path.exists(path):
                print(f"  MISSING: {path}")
                continue

            if scene not in cache:
                print(f"loading {scene} ...")
                with tifffile.TiffFile(path) as t:
                    cache[scene] = t.pages[0].asarray()

            a = cache[scene]
            if y + h > a.shape[0] or x + w > a.shape[1]:
                print(f"  OUT OF BOUNDS: {scene} {label} "
                      f"(image is {a.shape[1]}x{a.shape[0]})")
                continue

            crop = a[y:y+h, x:x+w]
            out = os.path.join(CROPS, f"{scene}__{label}.tif")
            tifffile.imwrite(out, crop)

            v = crop[crop > 0].astype(float)
            print(f"  {label:12s} {w}x{h}  mean={v.mean():6.1f}  "
                  f"CV={v.std()/v.mean():.3f}  -> {out}")
            written += 1

    print(f"\n{written} crops written to {CROPS}/")
    print("These are plain uint8 TIFFs. MATLAB: A = double(imread('...'));")


# ----------------------------------------------------------------------
if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "survey"
    if mode == "survey":
        survey()
    elif mode == "extract":
        extract()
    else:
        print(__doc__)
