#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

from PIL import Image
from scipy.ndimage import median_filter
from scipy.signal import savgol_filter


# ============================================================
# Configuration
# ============================================================

IMAGE = "input2.png"
OUTPUT = "digitized2.csv"

# Plot region in pixel coordinates.
# Top is clipped below the "DG:18 / PRF:3.9k / Filter:61" header text.
XMIN_PX = 15
XMAX_PX = 1500
YMIN_PX = 50
YMAX_PX = 464

# ------------------------------------------------------------
# Axis calibration
#
# Tick locations measured from the tick marks of the image:
#
# X (bottom axis, minor ticks every 56.15 px = 0.2 units):
#   pixel 103  -> -5
#   pixel 1507 ->  0
#
# Y (right axis, minor ticks every 41.6 px = 12.5 cm/s):
#   pixel 45.5  -> 100 cm/s
#   pixel 211.5 ->  50 cm/s
#   pixel 378.5 ->   0 cm/s   (matches the white baseline at row 377)
# ------------------------------------------------------------

X_PIX = np.array([103, 1507], dtype=float)
X_VAL = np.array([-5, 0], dtype=float)

Y_PIX = np.array([45.5, 211.5, 378.5], dtype=float)
Y_VAL = np.array([100, 50, 0], dtype=float)


# ============================================================
# Load image
# ============================================================

image = np.asarray(Image.open(IMAGE).convert("RGB"))

R = image[:, :, 0].astype(float)
G = image[:, :, 1].astype(float)
B = image[:, :, 2].astype(float)

gray = (R + G + B) / 3.0


# ============================================================
# Detect the grayscale Doppler spectrum
# ============================================================

# This trace is white/gray on a nearly black background, so brightness
# is the discriminator.  The teal caliper overlay (dotted line and the
# "A" marker) is rejected by requiring a roughly neutral colour.
BRIGHTNESS = 55

mask = (gray > BRIGHTNESS) & ((G - R) < 25) & ((B - R) < 25)

# Restrict detection to the actual plotting area.
mask[:YMIN_PX, :] = False
mask[YMAX_PX + 1:, :] = False
mask[:, :XMIN_PX] = False
mask[:, XMAX_PX + 1:] = False


# ============================================================
# Find the maximum-velocity envelope (signed)
# ============================================================

# Clinical convention: the digitized waveform is the *peak* (maximum)
# velocity envelope -- the outer edge of the spectrum -- and it is signed,
# so it follows the top edge during forward flow and the bottom edge while
# the flow is reversed.  It is never the bottom envelope alone, and never
# the middle of the band (that is the intensity-weighted mean velocity, a
# different quantity).
#
# Require a short vertical run of lit pixels so that isolated speckle
# outside the spectrum is not mistaken for the envelope.
RUN = 3

BASELINE_PX = 378.5

# A column counts as reverse flow when the spectrum extends this far below
# the zero line, for at least this many consecutive columns.  Both notches
# in this image are near-vertical zero crossings: a single pixel column
# holds the whole swing from about +50 to -13 cm/s, so simply taking the
# larger of the two extents picks the forward edge and misses the reversal.
# Detecting the reverse *window* first and then tracing its lower edge is
# what keeps the crossing intact.
REVERSE_PX = 12
REVERSE_MIN_COLS = 4

x_pixels = np.arange(XMIN_PX, XMAX_PX + 1)

top_y = np.full(len(x_pixels), np.nan)
bottom_y = np.full(len(x_pixels), np.nan)

kernel = np.ones(RUN, dtype=int)

for i, x in enumerate(x_pixels):

    run = np.convolve(mask[:, x].astype(int), kernel, mode="valid")

    ys = np.where(run == RUN)[0]

    if len(ys) == 0:
        continue

    top_y[i] = ys.min()
    bottom_y[i] = ys.max() + RUN - 1


# ------------------------------------------------------------
# Decide which side of the baseline each column belongs to
# ------------------------------------------------------------

reverse = np.nan_to_num(bottom_y - BASELINE_PX) > REVERSE_PX

# Drop reverse "windows" too short to be anything but speckle.
edges = np.diff(np.concatenate(([0], reverse.astype(int), [0])))
starts = np.where(edges == 1)[0]
stops = np.where(edges == -1)[0]

for start, stop in zip(starts, stops):
    if stop - start < REVERSE_MIN_COLS:
        reverse[start:stop] = False

envelope_y = np.where(reverse, bottom_y, top_y)

n_reverse = int(reverse.sum())
print(f"{n_reverse} of {len(x_pixels)} columns traced as reverse flow")


# ============================================================
# Fill missing points
# ============================================================

valid = np.isfinite(envelope_y)

if valid.sum() < 10:
    raise RuntimeError(
        "Could not detect enough bright pixels. "
        "Try lowering BRIGHTNESS."
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
envelope_y = median_filter(envelope_y, size=5)

# Savitzky-Golay gives a smooth waveform while preserving peaks.
envelope_y = savgol_filter(
    envelope_y,
    window_length=15,
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
plt.savefig("digitization_check2.png", dpi=200)


# ============================================================
# Plot physical data
# ============================================================

plt.figure(figsize=(12, 4))

plt.plot(x_data, velocity)

plt.xlabel("x")
plt.ylabel("Velocity [cm/s]")
plt.grid(True)

plt.tight_layout()

plt.show()
plt.savefig("digitized_waveform2.png", dpi=200)
