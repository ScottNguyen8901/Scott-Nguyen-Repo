import os
import numpy as np
import matplotlib.pyplot as plt

os.makedirs("plots", exist_ok=True)

# --------------------------------------------------
# Simple test image
# --------------------------------------------------
N = 256
img = np.zeros((N, N))
img[40:120, 40:120] = 1.0
img[150:220, 150:230] = 0.7
img[60:200:8, :] += 0.2
img = np.clip(img, 0, 1)

# --------------------------------------------------
# Fourier transform
# --------------------------------------------------
F = np.fft.fftshift(np.fft.fft2(img))

# Frequency grid
u = np.arange(N) - N // 2
v = np.arange(N) - N // 2
U, V = np.meshgrid(u, v)

# --------------------------------------------------
# 1) Rectangular low-pass
# --------------------------------------------------
rect = np.zeros((N, N))
rect[N//2-20:N//2+20, N//2-30:N//2+30] = 1

F_rect = F * rect
img_rect = np.real(np.fft.ifft2(np.fft.ifftshift(F_rect)))

plt.figure(figsize=(8, 3))
plt.subplot(1, 3, 1)
plt.imshow(img, cmap="gray")
plt.title("Original")
plt.axis("off")

plt.subplot(1, 3, 2)
plt.imshow(rect, cmap="gray")
plt.title("Rect LPF")
plt.axis("off")

plt.subplot(1, 3, 3)
plt.imshow(img_rect, cmap="gray")
plt.title("Output")
plt.axis("off")

plt.tight_layout()
plt.savefig("plots/problem_3_rectangular.png", dpi=300)
plt.close()

# --------------------------------------------------
# 2) Disk low-pass
# --------------------------------------------------
disk = ((U**2 + V**2) <= 25**2).astype(float)

F_disk = F * disk
img_disk = np.real(np.fft.ifft2(np.fft.ifftshift(F_disk)))

plt.figure(figsize=(8, 3))
plt.subplot(1, 3, 1)
plt.imshow(img, cmap="gray")
plt.title("Original")
plt.axis("off")

plt.subplot(1, 3, 2)
plt.imshow(disk, cmap="gray")
plt.title("Disk LPF")
plt.axis("off")

plt.subplot(1, 3, 3)
plt.imshow(img_disk, cmap="gray")
plt.title("Output")
plt.axis("off")

plt.tight_layout()
plt.savefig("plots/problem_3_disk.png", dpi=300)
plt.close()

# --------------------------------------------------
# 3) Circular support and rotation
# --------------------------------------------------
img_rot = np.rot90(img)

F_rot = np.fft.fftshift(np.fft.fft2(img_rot))
F_rot_disk = F_rot * disk
img_rot_disk = np.real(np.fft.ifft2(np.fft.ifftshift(F_rot_disk)))

# Rotate back to compare
img_rot_disk_back = np.rot90(img_rot_disk, -1)

plt.figure(figsize=(9, 3))
plt.subplot(1, 3, 1)
plt.imshow(img_disk, cmap="gray")
plt.title("Disk Filtered")
plt.axis("off")

plt.subplot(1, 3, 2)
plt.imshow(img_rot_disk_back, cmap="gray")
plt.title("Rotate Back")
plt.axis("off")

plt.subplot(1, 3, 3)
plt.imshow(np.abs(img_disk - img_rot_disk_back), cmap="gray")
plt.title("Difference")
plt.axis("off")

plt.tight_layout()
plt.savefig("plots/problem_3_rotation.png", dpi=300)
plt.close()

print("Saved:")
print("plots/problem_3_rectangular.png")
print("plots/problem_3_disk.png")
print("plots/problem_3_rotation.png")