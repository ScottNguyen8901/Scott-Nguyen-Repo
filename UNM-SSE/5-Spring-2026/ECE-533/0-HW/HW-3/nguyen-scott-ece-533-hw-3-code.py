
'''
Code to run main script
Usage:
    nguyen-scott-ece-533-hw-3-code.py 
'''
import os
import math
import time
import numpy as np
import cv2 as cv

from gabor_threads import build_filters, process, process_threaded
from dft import shift_dft


out_dir = "gabor_dft_outputs"
os.makedirs(out_dir, exist_ok=True)

size = 256
impulse = np.zeros((size, size), dtype=np.float64)
impulse[size // 2, size // 2] = 1.0

filters = build_filters()

for i, kern in enumerate(filters):
    filtered = cv.filter2D(impulse, cv.CV_64F, kern)

    h, w = filtered.shape[:2]
    dft_M = cv.getOptimalDFTSize(w)
    dft_N = cv.getOptimalDFTSize(h)

    dft_A = np.zeros((dft_N, dft_M, 2), dtype=np.float64)
    dft_A[:h, :w, 0] = filtered

    cv.dft(dft_A, dst=dft_A, nonzeroRows=h)

    re, im = cv.split(dft_A)
    magnitude = cv.sqrt(re ** 2.0 + im ** 2.0)
    log_spectrum = cv.log(1.0 + magnitude)

    shift_dft(log_spectrum, log_spectrum)
    cv.normalize(log_spectrum, log_spectrum, 0.0, 1.0, cv.NORM_MINMAX)

    save_path = os.path.join(out_dir, f"gabor_dft_{i:02d}.png")
    cv.imwrite(save_path, (log_spectrum * 255).astype(np.uint8))

print(f"Saved {len(filters)} DFT magnitude images to: {out_dir}")
cv.destroyAllWindows()


img = cv.imread("baboon.jpg")
if img is None:
    raise ValueError("Failed to load baboon.jpg (put it in the same folder).")

n = cv.getNumberOfCPUs()
print(f"Number of CPUs: {n}")

test_threads = [
    max(1, int(math.floor(0.5 * n))),
    max(1, int(n)),
    max(1, int(2 * n)),
]


start = time.time()
_ = process(img, filters)
T1 = time.time() - start
print(f"Single-thread time (16 filters): {T1:.4f} s")

results = {}
for k in test_threads:
    start = time.time()
    _ = process_threaded(img, filters, threadn=k)
    Tk = time.time() - start
    Sk = T1 / Tk
    results[k] = (Tk, Sk)
    print(f"Threads: {k:2d} | Time: {Tk:.4f} s | Speedup: {Sk:.3f}")

best_k = min(results, key=lambda kk: results[kk][0])
print(f"Best pool size (16 filters): {best_k} threads")


filters2 = filters + filters

start = time.time()
_ = process(img, filters2)
T1_2 = time.time() - start
print(f"\nSingle-thread time (32 filters): {T1_2:.4f} s")

results2 = {}
for k in test_threads:
    start = time.time()
    _ = process_threaded(img, filters2, threadn=k)
    Tk = time.time() - start
    Sk = T1_2 / Tk
    results2[k] = (Tk, Sk)
    print(f"Threads: {k:2d} | Time: {Tk:.4f} s | Speedup: {Sk:.3f}")

best_k2 = min(results2, key=lambda kk: results2[kk][0])
print(f"Best pool size (32 filters): {best_k2} threads")