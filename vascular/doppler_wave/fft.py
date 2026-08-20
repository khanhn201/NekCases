import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, detrend


radread = 0.004       # m
vis = 0.0035          # Pa s
rho = 1050.0          # kg/m^3
qmean = 0.250         # L/min


NMODES = 20


df = pd.read_csv('./digitized2.csv')
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





# ============================================================
# Split the record into individual cardiac cycles
#
# The six beats in this trace are NOT interchangeable: the R-R interval
# ranges from 0.78 to 0.98 s, and because systole stays roughly fixed
# while diastole stretches, the early-diastolic reverse notch lands
# anywhere between phase 0.10 and 0.22 once beats are normalized
# peak-to-peak.  The notch is only ~5% of a cycle wide, so averaging the
# beats -- which is what a periodic least-squares fit over the whole
# record does -- smears it away and leaves a waveform that never goes
# negative.  Fitting one representative beat keeps the triphasic shape.
# ============================================================

dx = np.mean(np.diff(x))

peaks, _ = find_peaks(
    flowrate,
    height=0.5 * (flowrate.max() + np.median(flowrate)),
    distance=int(0.6 / dx),
)

if len(peaks) < 2:
    raise ValueError("Could not segment the record into cardiac cycles.")

periods = np.diff(x[peaks])

print(f"\nFound {len(periods)} cycles, periods (s): "
      + " ".join(f"{p:.3f}" for p in periods))

# Resample every beat onto a common phase grid so they can be compared.
NPHASE = 200

phase = np.arange(NPHASE) / NPHASE

cycles = np.array([
    np.interp(
        phase,
        np.linspace(0.0, 1.0, stop - start, endpoint=False),
        flowrate[start:stop],
    )
    for start, stop in zip(peaks[:-1], peaks[1:])
])

print(" beat   period      min        max")
for i, (p, cyc) in enumerate(zip(periods, cycles)):
    print(f"{i+1:5d}  {p:7.3f}  {cyc.min(): .4f}  {cyc.max(): .4f}")


# ------------------------------------------------------------
# Pick the representative beat
#
# Among the beats that actually show flow reversal, take the one whose
# period is closest to the median -- i.e. a typical-length beat that is
# genuinely triphasic.  Set BEAT to a 1-based index to override.
# ------------------------------------------------------------

BEAT = None

triphasic = np.where(cycles.min(axis=1) < 0.0)[0]

if BEAT is not None:
    chosen = BEAT - 1
elif len(triphasic) > 0:
    chosen = triphasic[
        np.argmin(np.abs(periods[triphasic] - np.median(periods)))
    ]
else:
    chosen = int(np.argmin(np.abs(periods - np.median(periods))))

cycle = cycles[chosen]
wavelength = periods[chosen]

print(f"\nUsing beat {chosen+1} (period {wavelength:.4f} s, "
      f"min {cycle.min():.4f} L/min)")

# Rotate the cycle so that it starts at end-diastole (the foot of the
# systolic upstroke) rather than at the systolic peak.
foot = int(np.argmin(cycle[int(0.6 * NPHASE):]) + 0.6 * NPHASE)

cycle = np.roll(cycle, -foot)


# ============================================================
# Renormalize the chosen cycle to the target mean flow
# ============================================================

cycle_mean = np.mean(cycle)

cycle *= qmean / cycle_mean

n = NPHASE

print(f"Cycle mean flow after renormalization: {np.mean(cycle):.8f} L/min")


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
    cycle,
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

cycle_reconstructed = M @ coef

# ============================================================
# Frequency normalization
#
# The x axis of the Doppler strip is in seconds, so the measured heart
# rate comes straight from the beat period.  The chosen cycle is then
# time-rescaled to TARGET_BPM, i.e. the beat is stretched to a period of
# 60/TARGET_BPM seconds before the Fourier coefficients are interpreted
# as a frequency.
#
# The coefficients themselves are unchanged by this: the fit is done on
# a normalized phase grid, so rescaling time only reassigns which
# physical frequency phase = 1 corresponds to.  What it does change is
# the Womersley number and the stroke volume per beat.
#
# Set TARGET_BPM to None to keep the measured rate.
# ============================================================

TARGET_BPM = 60.0

measured_bpm = 60.0 / wavelength

if TARGET_BPM is None:
    freqbpm = measured_bpm
else:
    freqbpm = TARGET_BPM

freq_hz = freqbpm / 60.0

period = 1.0 / freq_hz

time_scale = period / wavelength

print(f"\nMeasured heart rate: {measured_bpm:.2f} bpm "
      f"(period {wavelength:.4f} s)")
print(f"Normalized to:       {freqbpm:.2f} bpm "
      f"(period {period:.4f} s), time scaled by {time_scale:.6f}")

# Stroke volume is mean flow times the (now normalized) beat period.
stroke_ml = qmean / 60.0 * period * 1000.0

print(f"Stroke volume at {freqbpm:.1f} bpm and {qmean:.3f} L/min: "
      f"{stroke_ml:.3f} mL/beat "
      f"(measured beat carried {qmean/60.0*wavelength*1000.0:.3f} mL)")

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
# Plot the chosen cycle and its reconstruction
# ============================================================

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 9))

for i, cyc in enumerate(cycles):
    ax1.plot(phase, cyc * qmean / cycle_mean, lw=0.8, alpha=0.5,
             label=f"beat {i+1}")

ax1.axhline(0.0, color="0.4", lw=1)
ax1.set_xlabel("Normalized phase (peak to peak)")
ax1.set_ylabel("Flow rate (L/min)")
ax1.set_title("All cycles — the reverse notch drifts between phase 0.10 "
              "and 0.22,\nwhich is why averaging them removes it")
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=8, ncol=3)

ax2.plot(phase, cycle, ".", ms=4, label=f"beat {chosen+1} (digitized)")

ax2.plot(phase, cycle_reconstructed, "-", lw=2,
         label=f"{NMODES}-mode Fourier reconstruction")

ax2.axhline(0.0, color="0.4", lw=1)
ax2.axhline(qmean, linestyle="--", label=f"Mean = {qmean:.3f} L/min")

ax2.set_xlabel(f"Phase (end-diastole to end-diastole, "
               f"1 phase unit = {period:.3f} s)")
ax2.set_ylabel("Flow rate (L/min)")
ax2.set_title(f"Representative cycle — {NMODES} Fourier modes\n"
              f"measured {measured_bpm:.1f} bpm, normalized to "
              f"{freqbpm:.1f} bpm")
ax2.grid(True, alpha=0.3)
ax2.legend()

plt.tight_layout()
plt.savefig("fft_check.png", dpi=150)
plt.show()
