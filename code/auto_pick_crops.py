#!/usr/bin/env python3
"""
auto_pick_crops.py — find real-data crops automatically

    python code/auto_pick_crops.py propose
        Scans every scene for 512 x 512 windows of three kinds, ranks them,
        and writes contact sheets to results\\figures\\candidates\\ :
          boundary : one vertical step between two regions   (Fig 6a analogue)
          pair     : a narrow vertical strip, two close edges (Fig 6e analogue)
          urban    : dense bright scatterers                  (Fig 6b analogue)

    python code/auto_pick_crops.py accept boundary=3 pair=1 urban=2
        Appends the chosen candidates to crops.csv with labels real_boundary,
        real_pair, real_urban. Then run extract_crops.py extract as usual.

Run from the repository root.

WHY VERTICAL: ROSED runs in the x-direction only, so it detects edges that
run top to bottom. Every score here rewards vertical structure and penalises
horizontal structure.

WHY YOU STILL LOOK: land and water cannot be separated reliably from simple
statistics (wind-roughened water can be brighter than land). The script ranks;
you reject any candidate with water on either side.

Selection only. Scoring uses a downsampled pyramid level and a log transform,
which is fine for FINDING regions. The detector itself always runs on the
untouched full-resolution crop.
"""
import os, sys, glob, json, csv
import numpy as np
import tifffile
from PIL import Image, ImageDraw

RAW, CROPFILE = "data/raw", "data/crops.csv"
OUT = "results/figures/candidates"
CAND = os.path.join(OUT, "candidates.json")
EXCLUDE = ["2025-06-22-23-57-52_UMBRA-10"]      # maritime, excluded (F8)
CROP = 512
TARGET_SCALE = 4                                 # read the ~1/4 pyramid level
TOP_PER_TYPE, MAX_PER_SCENE = 8, 3


def read_level(path):
    with tifffile.TiffFile(path) as t:
        H, W = t.pages[0].shape[:2]
        best = min(t.pages, key=lambda p: abs(H / p.shape[0] - TARGET_SCALE))
        a = best.asarray().astype(np.float64)
    return a, H, W, H / a.shape[0], W / a.shape[1]


def best_step(p, lo, hi):
    """Largest |mean(left) - mean(right)| over split positions lo..hi."""
    c = np.cumsum(p); n = len(p); best = (0.0, lo)
    for k in range(lo, hi):
        d = abs(c[k - 1] / k - (c[-1] - c[k - 1]) / (n - k))
        if d > best[0]: best = (d, k)
    return best


def band_score(p, w, flank):
    """Best narrow-band contrast: band of width w vs equal flanks either side."""
    n = len(p); best = (0.0, 0, 0.0)
    for x in range(flank + 24, n - w - flank - 24):
        L = p[x - flank:x].mean(); R = p[x + w:x + w + flank].mean()
        B = p[x:x + w].mean(); con = abs(B - (L + R) / 2)
        if abs(L - R) > 0.3 * con: continue           # flanks must match
        if max(p[x - flank:x].std(), p[x + w:x + w + flank].std()) > 0.5 * con:
            continue                                   # flanks must be quiet
        if con > best[0]: best = (con, x, abs(L - R))
    return best


def calib_boxes(scale_x, scale_y, scene):
    """Calibration regions from crops.csv, in preview coordinates."""
    boxes = []
    if not os.path.exists(CROPFILE): return boxes
    with open(CROPFILE) as fh:
        for r in csv.DictReader(fh):
            if r["scene"].strip() != scene: continue
            if not r["label"].strip().startswith(("clean", "calib")): continue
            x, y, w, h = (int(r[k]) for k in ("x", "y", "width", "height"))
            boxes.append((x / scale_x, y / scale_y, (x + w) / scale_x, (y + h) / scale_y))
    return boxes


def overlaps(b, boxes):
    x0, y0, x1, y1 = b
    return any(not (x1 <= a0 or a1 <= x0 or y1 <= b0 or b1 <= y0) for a0, b0, a1, b1 in boxes)


def scan_scene(path):
    name = os.path.basename(path)[:-4]
    a, H, W, sy, sx = read_level(path)
    L = np.log10(a + 1.0)
    win = int(round(CROP / sx)); stride = max(8, win // 4)
    p98 = np.percentile(L, 98)
    cal = calib_boxes(sx, sy, name)
    out = {"boundary": [], "pair": [], "urban": []}
    for py in range(0, L.shape[0] - win, stride):
        for px in range(0, L.shape[1] - win, stride):
            box = (px, py, px + win, py + win)
            if overlaps(box, cal): continue
            w = L[py:py + win, px:px + win]
            if (a[py:py + win, px:px + win] <= 0).mean() > 0.02: continue   # no-data
            col = w.mean(0); row = w.mean(1)
            colT = w[:win // 2].mean(0); colB = w[win // 2:].mean(0)
            lo, hi = win // 5, win - win // 5
            # --- boundary: one clean vertical step
            s, k = best_step(col, lo, hi)
            sT, kT = best_step(colT, lo, hi); sB, kB = best_step(colB, lo, hi)
            sr, _ = best_step(row, lo, hi)
            if s > 0 and abs(kT - kB) <= 3 and sr < 0.5 * s:
                resid = max(col[:k].std(), col[k:].std())
                if resid < 0.35 * s:
                    out["boundary"].append(dict(score=float(s), px=px, py=py, at=int(k)))
            # --- pair: narrow vertical strip with matching flanks
            bestp = (0.0, 0, 0.0, 0)
            for bw in range(3, 9):                     # 12-32 full-res px at scale 4
                c, x, d = band_score(col, bw, 12)
                if c > bestp[0]: bestp = (c, x, d, bw)
            if bestp[0] > 0:
                cT, xT, _ = band_score(colT, bestp[3], 12)
                cB, xB, _ = band_score(colB, bestp[3], 12)
                if cT > 0 and cB > 0 and abs(xT - xB) <= 2:
                    out["pair"].append(dict(score=float(bestp[0]), px=px, py=py,
                                            at=int(bestp[1]), width=int(bestp[3])))
            # --- urban: many bright scatterers, vertically organised
            bright = (w > p98).mean()
            if bright > 0.08 and col.std() > 1.2 * row.std():
                out["urban"].append(dict(score=float(bright), px=px, py=py))
    for t in out:
        for c in out[t]:
            c.update(scene=name, sx=sx, sy=sy, W=W, H=H, win=win)
    return name, a, out


def nms(cands, win):
    cands = sorted(cands, key=lambda c: -c["score"]); keep, per = [], {}
    for c in cands:
        if per.get(c["scene"], 0) >= MAX_PER_SCENE: continue
        if any(k["scene"] == c["scene"] and abs(k["px"] - c["px"]) < win and
               abs(k["py"] - c["py"]) < win for k in keep): continue
        keep.append(c); per[c["scene"]] = per.get(c["scene"], 0) + 1
        if len(keep) >= TOP_PER_TYPE: break
    return keep


def stretch(x):
    lo, hi = np.percentile(x, [1, 99])
    return (np.clip((x - lo) / max(hi - lo, 1e-9), 0, 1) * 255).astype(np.uint8)


def propose():
    os.makedirs(OUT, exist_ok=True)
    files = [f for f in sorted(glob.glob(os.path.join(RAW, "*_GEC.tif")))
             if not any(e in f for e in EXCLUDE)]
    if not files: sys.exit(f"No scenes in {RAW}")
    allc = {"boundary": [], "pair": [], "urban": []}; imgs = {}
    for f in files:
        print("scanning", os.path.basename(f), flush=True)
        name, a, out = scan_scene(f); imgs[name] = a
        for t in out: allc[t] += out[t]
    chosen = {}
    for t in allc:
        win = allc[t][0]["win"] if allc[t] else 128
        chosen[t] = nms(allc[t], win)
        for i, c in enumerate(chosen[t], 1):
            c["id"] = i
            c["x0"] = int(min(max(round(c["px"] * c["sx"]), 0), c["W"] - CROP))
            c["y0"] = int(min(max(round(c["py"] * c["sy"]), 0), c["H"] - CROP))
        # contact sheet: each candidate at 3x, with its id and scene
        tiles = []
        for c in chosen[t]:
            w = c["win"]; a = imgs[c["scene"]][c["py"]:c["py"] + w, c["px"]:c["px"] + w]
            im = Image.fromarray(stretch(a)).resize((w * 3, w * 3), Image.NEAREST).convert("RGB")
            d = ImageDraw.Draw(im)
            d.rectangle([0, 0, 150, 30], fill=(0, 0, 0))
            d.text((6, 8), f"{t} {c['id']}", fill=(255, 255, 0))
            d.text((6, w * 3 - 16), c["scene"][:19], fill=(255, 255, 0))
            tiles.append(im)
        if tiles:
            tw, th = tiles[0].size; cols = 4; rows = (len(tiles) + cols - 1) // cols
            sheet = Image.new("RGB", (cols * tw, rows * th), (40, 40, 40))
            for i, im in enumerate(tiles):
                sheet.paste(im, ((i % cols) * tw, (i // cols) * th))
            sheet.save(os.path.join(OUT, f"sheet_{t}.png"))
        print(f"\n{t}: {len(chosen[t])} candidates -> {OUT}/sheet_{t}.png")
        for c in chosen[t]:
            extra = f"  strip width ~{c['width']*4} px" if t == "pair" else ""
            print(f"  {c['id']}  {c['scene']}  x={c['x0']} y={c['y0']}  score {c['score']:.3f}{extra}")
    for t in chosen:
        for c in chosen[t]: c.pop("W", None)
    json.dump(chosen, open(CAND, "w"), indent=1)
    print("\nOpen the three sheet_*.png files. Reject any tile with water on either")
    print("side or edges that do not run top to bottom. Then, for example:")
    print("  python code/auto_pick_crops.py accept boundary=1 pair=2 urban=1")


def accept(args):
    if not os.path.exists(CAND): sys.exit("Run 'propose' first.")
    chosen = json.load(open(CAND)); lines = []
    for a in args:
        t, i = a.split("="); c = [c for c in chosen[t] if c["id"] == int(i)][0]
        lines.append(f"{c['scene']},real_{t},{c['x0']},{c['y0']},{CROP},{CROP}")
    with open(CROPFILE) as fh: text = fh.read()
    keep = [ln for ln in text.splitlines() if ln.strip() and
            not any(f",real_{a.split('=')[0]}," in ln for a in args)]
    with open(CROPFILE, "w") as fh: fh.write("\n".join(keep + lines) + "\n")
    print("Added to crops.csv (replacing any earlier rows with the same labels):")
    for ln in lines: print("  " + ln)
    print("\nNext:  python code/extract_crops.py extract")


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in ("propose", "accept"): print(__doc__)
    elif sys.argv[1] == "propose": propose()
    else: accept(sys.argv[2:])
