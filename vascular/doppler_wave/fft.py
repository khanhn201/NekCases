import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, detrend


radread = 0.004       # m
vis = 0.0035          # Pa s
rho = 1050.0          # kg/m^3
qmean = 0.250         # L/min


NMODES = 20


df = pd.read_csv('./digitized.csv')
x = df.x.to_numpy()
velocity = df["velocity_cm_s"].to_numpy(dtype=float)


# ============================================================
# Convert velocity [cm/s] -> flow rate [L/min]
#
# radius = 0.004 m = 0.4 cm
#
# Q [cm^3/s] = velocity [cm/s] * area [cm^2]
# 1 L = 1000 cm^3
# ============================================================
radius_cm = radread * 100.0
area_cm2 = np.pi * radius_cm**2
flowrate = velocity * area_cm2       # cm^3/s
flowrate = flowrate * 60.0 / 1000.0  # L/min
print(f"Original mean flow: {np.mean(flowrate):.8f} L/min")

original_mean = np.mean(flowrate)
if original_mean == 0:
    raise ValueError("Mean flow rate is zero; cannot normalize.")
scale = qmean / original_mean
flowrate *= scale
print(f"Normalization factor: {scale:.8f}")
print(f"Normalized mean flow: {np.mean(flowrate):.8f} L/min")





# Frequency is detected from the dominant spectral peak of the detrended signal.
# Since x is spatial rather than time, normalize the dominant spatial period to 1.
dx = np.mean(np.diff(x))
n = len(flowrate)

q_detrended = flowrate - np.mean(flowrate)

freq = np.fft.rfftfreq(n, d=dx)
fft = np.fft.rfft(q_detrended)

power = np.abs(fft)

# Remove DC component
power[0] = 0.0

k0 = np.argmax(power)

f0 = abs(freq[k0])

if f0 <= 0:
    raise ValueError("Could not detect a nonzero fundamental frequency.")

wavelength = 1.0 / f0

print(f"Detected fundamental frequency: {f0:.10e} cycles/x")
print(f"Detected period/wavelength:      {wavelength:.10e} x-units")



phase = (x - x[0]) / wavelength

# ============================================================
# 20-mode Fourier least-squares decomposition
#
# Q(t) = Qmean
#      + sum [ A_k cos(2*pi*k*phase)
#             +B_k sin(2*pi*k*phase) ]
# ============================================================

M = np.ones((n, 2 * NMODES + 1))

for k in range(1, NMODES + 1):
    M[:, 2*k - 1] = np.cos(2.0 * np.pi * k * phase)
    M[:, 2*k]     = np.sin(2.0 * np.pi * k * phase)

coef, residuals, rank, singular_values = np.linalg.lstsq(
    M,
    flowrate,
    rcond=None
)

# Constant component
fourier_mean = coef[0]

# Fourier coefficients
A = coef[1::2]
B = coef[2::2]

print(f"Fourier mean: {fourier_mean:.12e} L/min")

# ============================================================
# Reconstruct 20-mode waveform
# ============================================================

flowrate_reconstructed = M @ coef

# ============================================================
# Frequency for Womersley calculation
#
# The cl.fft routine expects freqbpm.
#
# We normalized the independent variable to one fundamental
# period, so the physical frequency used here is 1 Hz.
# ============================================================

freq_hz = 1.0
freqbpm = 60.0

omega = 2.0 * np.pi * freq_hz

womp = radread * np.sqrt(
    omega * rho / vis
)

# ============================================================
# Write cl.fft
# ============================================================

jmax = 0
imax = 0

start = 0.0
end = 1.0

itype = 1

with open('cl.fft', "w") as f:

    # Header
    f.write(
        f"{jmax:d} {imax:d} "
        f"{start:.16e} {end:.16e} "
        f"{vis:.16e} {rho:.16e} "
        f"{radread:.16e} {freqbpm:.16e} "
        f"{womp:.16e}\n"
    )

    # Fourier information
    f.write(
        f"{itype:d} {NMODES:d} {qmean:.16e}\n"
    )

    # A_k B_k
    for k in range(NMODES):
        f.write(
            f"{A[k]:.16e} {B[k]:.16e}\n"
        )

print(f"\nWrote cl.fft")

# ============================================================
# Print coefficients
# ============================================================

print("\nFourier coefficients:")
print(" k             A_k                  B_k")

for k in range(NMODES):
    print(
        f"{k+1:2d}  "
        f"{A[k]: .12e}  "
        f"{B[k]: .12e}"
    )

# ============================================================
# Plot original normalized waveform and reconstruction
# ============================================================

order = np.argsort(phase)

plt.figure(figsize=(10, 5))

plt.plot(
    phase,
    flowrate,
    ".",
    ms=3,
    label="Normalized flow rate"
)

plt.plot(
    phase[order],
    flowrate_reconstructed[order],
    "-",
    lw=2,
    label="20-mode Fourier reconstruction"
)

plt.axhline(
    qmean,
    linestyle="--",
    label=f"Mean = {qmean:.3f} L/min"
)

plt.xlabel("Normalized phase")
plt.ylabel("Flow rate (L/min)")
plt.title(
    f"Normalized flow waveform — 20 Fourier modes\n"
    f"Fundamental frequency = {freq_hz:.1f} Hz"
)

plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()
