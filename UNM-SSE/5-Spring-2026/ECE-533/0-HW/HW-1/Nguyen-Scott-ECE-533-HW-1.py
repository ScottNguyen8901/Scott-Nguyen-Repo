"""
=====================================================
ECE 533 - Homework 1
Scott Nguyen
Problems 1–6
All outputs saved to ./plots
=====================================================
"""

import os
import numpy as np
import cv2
import matplotlib.pyplot as plt
import torch
import torch.nn.functional as F
from AOLME import *
from IPython.display import HTML


# =================================================
# Functions
# =================================================
def ensure_plot_dir():
    base = os.path.dirname(__file__)
    p = os.path.join(base, "plots")
    os.makedirs(p, exist_ok=True)
    return p


def save_gray(img, title, filename, vmin=0, vmax=255):
    plt.figure()
    plt.imshow(img, cmap='gray', vmin=vmin, vmax=vmax)
    plt.title(title)
    plt.axis('off')
    plt.savefig(os.path.join(plot_dir, filename), dpi=300, bbox_inches='tight')
    plt.close()


def save_rgb(img, title, filename):
    plt.figure()
    plt.imshow(img)
    plt.title(title)
    plt.axis('off')
    plt.savefig(os.path.join(plot_dir, filename), dpi=300, bbox_inches='tight')
    plt.close()


def log_transform(I, C):
    J = np.log1p(C * I.astype(np.float32))
    J = (J - J.min()) / (J.max() - J.min()) * 255
    return J.astype(np.uint8)


def relu(x):
    return np.maximum(0, x)


def sigmoid(x):
    return 1 / (1 + np.exp(-x))


def max_pool(I, k, stride=1):
    t = torch.from_numpy(I).float()[None, None, :, :]
    return F.max_pool2d(t, kernel_size=k, stride=stride).squeeze().numpy()


def erosion(I, k):
    t = torch.from_numpy(I).float()[None, None, :, :]
    return (-F.max_pool2d(-t, kernel_size=k, stride=1)).squeeze().numpy()


def imread_rgb_resize(path, h=320):
    bgr = cv2.imread(path)
    H, W = bgr.shape[:2]
    bgr = cv2.resize(bgr, (int(W * h / H), h))
    return cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB)


def threshold_mask(rgb, r0, r1, g0, g1, b0, b1):
    r, g, b = rgb[:,:,0], rgb[:,:,1], rgb[:,:,2]
    return ((r>=r0)&(r<=r1)&(g>=g0)&(g<=g1)&(b>=b0)&(b<=b1)).astype(np.uint8)


def apply_mask(rgb, m):
    out = rgb.copy()
    out[m==0] = 0
    return out


# =================================================
plot_dir = ensure_plot_dir()


# =================================================
# Problem 1
# =================================================
binary = np.indices((8,8)).sum(axis=0) % 2
save_gray(binary, "Problem 1(a) Binary", "binary.png", 0, 1)

row = np.linspace(0,255,8,dtype=np.uint8)
gray = np.tile(row,(8,1))
save_gray(gray, "Problem 1(b) Grayscale", "grayscale.png")

color = np.zeros((8,8,3),dtype=np.uint8)
color[:4,:4]=[255,0,0]
color[:4,4:]=[0,255,0]
color[4:,:4]=[0,0,255]
color[4:,4:]=[255,255,0]
save_rgb(color, "Problem 1(c) Color", "color.png")


# =================================================
# Problem 2
# =================================================
A = np.full((128,128),128,np.uint8)
A[48:80,48:80]=129
save_gray(A,"P2(A) Min contrast","P2A.png")

B = np.zeros((128,128),np.uint8)
B[32:96,32:96]=255
save_gray(B,"P2(B) Max contrast","P2B.png")

vals=[5,10,100,200,250]
C = np.zeros((200,400),np.uint8)
for i,v in enumerate(vals):
    C[40:160,30+i*80:90+i*80]=v
save_gray(C,"P2(C) Patches","P2C.png")

for c in [1,5,10,100]:
    save_gray(log_transform(C,c),f"P2(D) Log C={c}",f"P2D_{c}.png")


# =================================================
# Problem 3
# =================================================
img=np.zeros((150,300),np.uint8)
img[40:120,40:120]=5
img[40:120,180:260]=100
save_gray(img,"P3(a) Original","P3A.png")

rimg = relu(img-20)
rimg=(rimg/rimg.max()*255).astype(np.uint8)
save_gray(rimg,"P3(b) ReLU","P3B.png")

simg=sigmoid(0.1*(img-20))
simg=((simg-simg.min())/(simg.max()-simg.min())*255).astype(np.uint8)
save_gray(simg,"P3(c) Sigmoid","P3C.png")


# =================================================
# Problem 4
# =================================================
I=np.random.randint(0,256,(64,64)).astype(np.float32)
save_gray(I,"P4 Input","P4A.png")

pool=max_pool(I,4,4)
save_gray(pool,"P4 MaxPool","P4B.png")

ero=erosion(I,3)
save_gray(ero,"P4 Erosion","P4C.png")


# =================================================
# Problem 5
# =================================================
rgb=imread_rgb_resize("Bob.jpg")

cases={
"a_bear":(0,120,0,120,0,120),
"b_skin":(120,255,120,255,0,160),
"c_clothes":(0,120,0,120,120,255)
}

for name,t in cases.items():
    m=threshold_mask(rgb,*t)
    out=apply_mask(rgb,m)
    save_rgb(rgb,f"P5 {name} input",f"P5_{name}_input.png")
    save_gray(m*255,f"P5 {name} mask",f"P5_{name}_mask.png")
    save_rgb(out,f"P5 {name} output",f"P5_{name}_output.png")


# =================================================
# Problem 6
# =================================================
rows,cols=20,20
gray="%02x%02x%02x"%(120,120,120)
orange="%02x%02x%02x"%(255,165,0)
cyan="%02x%02x%02x"%(0,255,255)

frames=[]
r1,c1,dr1,dc1=2,2,1,1
r2,c2,dr2,dc2=14,12,-1,1

for _ in range(20):
    f=np.full((rows,cols),gray)
    f[r1:r1+2,c1:c1+2]=orange
    f[r2:r2+3,c2:c2+3]=cyan
    frames.append(f)

    r1+=dr1; c1+=dc1
    r2+=dr2; c2+=dc2

    if r1<=0 or r1>=rows-2: dr1*=-1
    if c1<=0 or c1>=cols-2: dc1*=-1
    if r2<=0 or r2>=rows-3: dr2*=-1
    if c2<=0 or c2>=cols-3: dc2*=-1

play_video = vid_show(frames, 10)
play_video.save(os.path.join(plot_dir, "problem6_video.mp4"))