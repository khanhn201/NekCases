#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

from PIL import Image
from scipy.ndimage import median_filter
from scipy.signal import savgol_filter


# ============================================================
# Configuration
# ============================================================

IMAGE = "input.png"
OUTPUT = "digitized.csv"

# Plot region in pixel coordinates.
# These values are tuned for the supplied image.
XMIN_PX = 8
XMAX_PX = 1120
YMIN_PX = 25
YMAX_PX = 330

# ------------------------------------------------------------
# Axis calibration
#
# Pixel locations measured from the supplied image:
#
# X:
#   pixel 77   -> -4
#   pixel 337  -> -3
#   pixel 597  -> -2
#   pixel 858  -> -1
#   pixel 1119 ->  0
#
# Y:
#   pixel 94   -> 100 cm/s
#   pixel 212  ->  50 cm/s
#   pixel 330  ->   0 cm/s
# ------------------------------------------------------------

X_PIX = np.array([77, 337, 597, 858, 1119], dtype=float)
X_VAL = np.array([-4, -3, -2, -1, 0], dtype=float)

Y_PIX = np.array([94, 212, 330], dtype=float)
Y_VAL = np.array([100, 50, 0], dtype=float)


# ============================================================
# Load image
# ============================================================

image = np.asarray(Image.open(IMAGE).convert("RGB"))

R = image[:, :, 0].astype(float)
G = image[:, :, 1].astype(float)
B = image[:, :, 2].astype(float)


# ============================================================
# Detect orange/red Doppler signal
# ============================================================

# The signal is orange/red, while the background is nearly black.
#
# These thresholds can be adjusted if necessary.
mask = (
    (R > 80)
    & (R > 1.15 * G)
    & (R > 1.25 * B)
    & ((R - G) > 20)
)

# Restrict detection to the actual plotting area.
mask[:YMIN_PX, :] = False
mask[YMAX_PX + 1:, :] = False
mask[:, :XMIN_PX] = False
mask[:, XMAX_PX + 1:] = False


# ============================================================
# Find upper envelope
# ============================================================

x_pixels = np.arange(XMIN_PX, XMAX_PX + 1)

envelope_y = np.full(len(x_pixels), np.nan)

for i, x in enumerate(x_pixels):

    ys = np.where(mask[:, x])[0]

    if len(ys) > 0:
        envelope_y[i] = ys.min()


# ============================================================
# Fill missing points
# ============================================================

valid = np.isfinite(envelope_y)

if valid.sum() < 10:
    raise RuntimeError(
        "Could not detect enough red/orange pixels. "
        "Try changing the color thresholds."
    )

envelope_y = np.interp(
    np.arange(len(envelope_y)),
    np.where(valid)[0],
    envelope_y[valid],
)


# ============================================================
# Smooth the envelope
# ============================================================

# Median filtering removes isolated pixel spikes.
envelope_y = median_filter(envelope_y, size=9)

# Savitzky-Golay gives a smooth waveform while preserving peaks.
envelope_y = savgol_filter(
    envelope_y,
    window_length=41,
    polyorder=3,
)


# ============================================================
# Pixel -> physical coordinate calibration
# ============================================================

# x = a*x_pixel + b
x_coeff = np.polyfit(X_PIX, X_VAL, 1)
x_data = np.polyval(x_coeff, x_pixels)

# velocity = a*y_pixel + b
y_coeff = np.polyfit(Y_PIX, Y_VAL, 1)
velocity = np.polyval(y_coeff, envelope_y)


# ============================================================
# Save data
# ============================================================

data = np.column_stack((x_data, velocity))

np.savetxt(
    OUTPUT,
    data,
    delimiter=",",
    header="x,velocity_cm_s",
    comments="",
)

print(f"Saved {len(data)} points to {OUTPUT}")


# ============================================================
# Diagnostic plot
# ============================================================

fig, ax = plt.subplots(figsize=(12, 5))

ax.imshow(image)

ax.plot(
    x_pixels,
    envelope_y,
    linewidth=1.5,
    label="Extracted envelope",
)

ax.set_xlim(XMIN_PX, XMAX_PX)
ax.set_ylim(YMAX_PX, YMIN_PX)
ax.set_xlabel("Image x [pixel]")
ax.set_ylabel("Image y [pixel]")
ax.legend()

plt.tight_layout()
plt.savefig("digitization_check.png", dpi=200)
plt.show()


# ============================================================
# Plot physical data
# ============================================================

plt.figure(figsize=(12, 4))

plt.plot(x_data, velocity)

plt.xlabel("x")
plt.ylabel("Velocity [cm/s]")
plt.grid(True)

plt.tight_layout()
plt.savefig("digitized_waveform.png", dpi=200)
plt.show()
