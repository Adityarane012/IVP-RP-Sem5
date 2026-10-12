import tifffile, numpy as np, glob, os
for f in sorted(glob.glob('data/crops/*__calib.tif')):
    a = tifffile.imread(f).astype(float)
    best = min(
        ((a[r:r+128, c:c+128].std()/max(a[r:r+128, c:c+128].mean(),1e-9), r, c)
         for r in range(0, 385, 32) for c in range(0, 385, 32)),
        key=lambda t: t[0])
    print(f'{os.path.basename(f)[:28]:30s} best CV={best[0]:.3f} at ({best[2]},{best[1]})')