#!/usr/bin/env python
'''
sample for discrete fourier transform (dft)

USAGE:
dft.py <image_file>
'''

from __future__ import print_function
import numpy as np
import cv2 as cv
import sys


def shift_dft(src, dst=None):
    '''
    Rearrange the quadrants of Fourier image so that the origin is at
    the image center. Swaps quadrant 1 with 3, and 2 with 4.
    src and dst arrays must be equal size & type
    '''
    if dst is None:
        dst = np.empty(src.shape, src.dtype)
    elif src.shape != dst.shape:
        raise ValueError("src and dst must have equal sizes")
    elif src.dtype != dst.dtype:
        raise TypeError("src and dst must have equal types")

    if src is dst:
        ret = np.empty(src.shape, src.dtype)
    else:
        ret = dst

    h, w = src.shape[:2]
    cx1 = cx2 = w // 2
    cy1 = cy2 = h // 2

    if w % 2 != 0:
        cx2 += 1
    if h % 2 != 0:
        cy2 += 1

    ret[h-cy1:, w-cx1:] = src[0:cy1, 0:cx1]
    ret[0:cy2, 0:cx2] = src[h-cy2:, w-cx2:]
    ret[0:cy2, w-cx2:] = src[h-cy2:, 0:cx2]
    ret[h-cy1:, 0:cx1] = src[0:cy1, w-cx1:]

    if src is dst:
        dst[:, :] = ret

    return dst


def main():
    if len(sys.argv) > 1:
        fname = sys.argv[1]
    else:
        fname = 'baboon.jpg'

    print("usage: python dft.py <image_file>")

    im = cv.imread(cv.samples.findFile(fname))
    im = cv.cvtColor(im, cv.COLOR_BGR2GRAY)

    h, w = im.shape[:2]
    realInput = im.astype(np.float64)

    dft_M = cv.getOptimalDFTSize(w)
    dft_N = cv.getOptimalDFTSize(h)

    dft_A = np.zeros((dft_N, dft_M, 2), dtype=np.float64)
    dft_A[:h, :w, 0] = realInput

    cv.dft(dft_A, dst=dft_A, nonzeroRows=h)

    cv.imshow("win", im)

    image_Re, image_Im = cv.split(dft_A)
    magnitude = cv.sqrt(image_Re ** 2.0 + image_Im ** 2.0)
    log_spectrum = cv.log(1.0 + magnitude)

    shift_dft(log_spectrum, log_spectrum)

    cv.normalize(log_spectrum, log_spectrum, 0.0, 1.0, cv.NORM_MINMAX)
    cv.imshow("magnitude", log_spectrum)

    cv.waitKey(0)
    print('Done')


if __name__ == '__main__':
    print(__doc__)
    main()
    cv.destroyAllWindows()