import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.special import jv


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

# Fourier coefficients of the ENVELOPE (centreline velocity) waveform.
# These are converted to flow-rate coefficients further down.
A_env = coef[1::2]
B_env = coef[2::2]

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
# Peak -> flow rate:  Womersley inversion
#
# The Doppler max envelope is the PEAK (centreline) velocity, not the
# cross-sectional mean, so the flow rate is not simply velocity * area.
# For steady flow the centreline/mean ratio is 2 (parabolic), but under
# pulsatile flow it is frequency dependent: the profile goes blunt as
# the Womersley number rises, so each harmonic has its own complex ratio
#
#     G_k = u_centreline,k / u_mean,k
#         = [1 - 1/J0(L)] / [1 - 2*J1(L)/(L*J0(L],   L = alpha_k * i^(3/2)
#
# with alpha_k = alpha_1 * sqrt(k).  G_0 = 2 (parabolic), and |G_k| falls
# to ~1.06 by k = 20.
#
# This matters because womer() in ab.usr rebuilds the Womersley profile
# from the flow coefficients.  Feeding it coefficients derived with a
# constant (plug or parabolic) factor means the centreline velocity the
# solver actually imposes does NOT reproduce the digitized envelope --
# it comes out ~1.4x less pulsatile, and the reverse-flow phase vanishes.
#
# Set WOMERSLEY_INVERSION = False to fall back to Q proportional to the
# envelope.
# ============================================================

WOMERSLEY_INVERSION = True

alpha1 = womp

# Complex harmonics of the envelope, in the convention used by ab.usr:
#   waveform(t) = mean + Re[ sum_k (A_k - i B_k) * exp(i k omega t) ]
v_hat = A_env - 1j * B_env

if WOMERSLEY_INVERSION:

    k_modes = np.arange(1, NMODES + 1)

    lam = alpha1 * np.sqrt(k_modes) * np.exp(0.75j * np.pi)

    G = (1.0 - 1.0 / jv(0, lam)) / (
        1.0 - 2.0 * jv(1, lam) / (lam * jv(0, lam))
    )

    # G_0 = 2 for the mean, so harmonics are rescaled by 2/G_k relative
    # to a mean that is held at qmean.
    q_hat = v_hat * 2.0 / G

else:
    G = np.full(NMODES, 2.0, dtype=complex)
    q_hat = v_hat.astype(complex)

A = q_hat.real
B = -q_hat.imag

print(f"\nWomersley inversion: {WOMERSLEY_INVERSION}  (alpha_1 = {alpha1:.4f})")

if WOMERSLEY_INVERSION:
    print("   k   alpha_k    |G_k|   arg(G_k)   harmonic gain 2/|G_k|")
    for k in (1, 2, 3, 5, 10, 20):
        g = G[k-1]
        print(f"{k:4d}  {alpha1*np.sqrt(k):7.2f}  {abs(g):7.3f}  "
              f"{np.degrees(np.angle(g)):7.2f} deg  {2.0/abs(g):8.3f}")

# ------------------------------------------------------------
# Waveforms for plotting / diagnostics
# ------------------------------------------------------------

flow_cycle = np.full(NPHASE, qmean)

for k in range(1, NMODES + 1):
    flow_cycle += (A[k-1] * np.cos(2.0 * np.pi * k * phase)
                   + B[k-1] * np.sin(2.0 * np.pi * k * phase))

# Round-trip: the centreline velocity womer() will build from these
# coefficients.  This is a check on the A/B <-> complex convention, not
# an independent validation -- it should match the envelope to within
# the 20-mode truncation.
centreline = np.full(NPHASE, qmean)

for k in range(1, NMODES + 1):
    c = (A[k-1] - 1j * B[k-1]) * G[k-1] / 2.0
    centreline += (c.real * np.cos(2.0 * np.pi * k * phase)
                   - c.imag * np.sin(2.0 * np.pi * k * phase))

print(f"\nround-trip centreline vs envelope: max |diff| = "
      f"{np.max(np.abs(centreline - cycle_reconstructed)):.3e} L/min-equiv")

fwd = flow_cycle[flow_cycle > 0].sum()
rev = abs(flow_cycle[flow_cycle < 0].sum())

print(f"\nFlow waveform written to cl.fft:")
print(f"   min {flow_cycle.min():.4f}  max {flow_cycle.max():.4f} L/min"
      f"   peak/mean {flow_cycle.max()/qmean:.2f}"
      f"   retrograde fraction {rev/fwd:.3f}")

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
# Axial velocity profiles the solver will actually impose
#
# womer() in ab.usr rebuilds u_z(r,t) from the coefficients just
# written to cl.fft:
#
#   u(rhat,t) = 2*Qmean/Area * (1 - rhat^2)
#             + Re{ sum_k  q_k/Area * PHI_k(rhat) * exp(i k omega t) }
#
#   PHI_k(rhat) = [1 - J0(lam_k*rhat)/J0(lam_k)]
#                 / [1 - 2*J1(lam_k)/(lam_k*J0(lam_k))]
#
# with lam_k = alpha_k * i^(3/2) and q_k = A_k - i B_k.  PHI_k(0) = G_k,
# so the centreline of these profiles is the envelope fitted above.  The
# Bessel kernel is used regardless of WOMERSLEY_INVERSION, because that
# is what womer() does with whatever coefficients it is handed.
# ============================================================

NR = 101

rhat = np.linspace(0.0, 1.0, NR)

area = np.pi * radread**2

# Q [L/min] -> [m^3/s]
lpm_to_m3s = 1.0 / (1000.0 * 60.0)

lam_k = alpha1 * np.sqrt(np.arange(1, NMODES + 1)) * np.exp(0.75j * np.pi)

# PHI[k-1, :] : radial shape of harmonic k
PHI = np.array([
    (1.0 - jv(0, lam * rhat) / jv(0, lam))
    / (1.0 - 2.0 * jv(1, lam) / (lam * jv(0, lam)))
    for lam in lam_k
])

q_k = (A - 1j * B) * lpm_to_m3s / area          # m/s, complex


def axial_velocity(t_phase):
    """Axial velocity profile u(rhat) [m/s] at normalized phase t_phase."""
    u = 2.0 * qmean * lpm_to_m3s / area * (1.0 - rhat**2)
    for k in range(1, NMODES + 1):
        u += np.real(q_k[k-1] * PHI[k-1] * np.exp(2j * np.pi * k * t_phase))
    return u


# ------------------------------------------------------------
# Phases to report
#
#   systole        peak forward flow
#   diastole       peak reverse flow (the early-diastolic notch)
#   late diastole  end-diastole, i.e. phase 0 -- the cycle was rolled so
#                  it starts at the foot of the systolic upstroke
# ------------------------------------------------------------

i_sys = int(np.argmax(flow_cycle))
i_dia = int(np.argmin(flow_cycle))
i_late = 0

phase_points = [
    ("diastole",      phase[i_dia],  flow_cycle[i_dia]),
    ("late_diastole", phase[i_late], flow_cycle[i_late]),
    ("systole",       phase[i_sys],  flow_cycle[i_sys]),
]

profiles = {name: axial_velocity(ph) for name, ph, _ in phase_points}

print("\nAxial velocity profiles imposed by womer() "
      f"(alpha_1 = {alpha1:.3f}, period {period:.3f} s):")
print("  phase name      phase    t [s]     Q [L/min]   "
      "u_centreline [m/s]   u_mean [m/s]")

for name, ph, q in phase_points:
    u = profiles[name]
    # cross-sectional mean = 2*int_0^1 u rhat drhat
    u_bar = 2.0 * np.trapezoid(u * rhat, rhat)
    print(f"  {name:<14s} {ph:6.3f}  {ph*period:7.4f}  {q: 10.4f}   "
          f"{u[0]: 14.6f}      {u_bar: 10.6f}")

# ------------------------------------------------------------
# Write the profiles
# ------------------------------------------------------------

with open("velocity_profiles.dat", "w") as f:
    f.write("# Axial velocity profiles rebuilt from cl.fft, as womer() in "
            "ab.usr imposes them\n")
    f.write(f"# radius R = {radread:.6e} m,  period = {period:.6f} s "
            f"({freqbpm:.1f} bpm),  alpha_1 = {alpha1:.6f}\n")
    f.write(f"# mean flow = {qmean:.6f} L/min,  "
            f"mean velocity = {qmean*lpm_to_m3s/area:.6f} m/s\n")
    for name, ph, q in phase_points:
        f.write(f"# {name}: phase = {ph:.4f}, t = {ph*period:.6f} s, "
                f"Q = {q:.6f} L/min\n")
    f.write("# columns: r/R  r[m]  "
            + "  ".join(f"u_{name}[m/s]" for name, _, _ in phase_points)
            + "\n")
    for i in range(NR):
        f.write(f"{rhat[i]:.16e} {rhat[i]*radread:.16e} "
                + " ".join(f"{profiles[name][i]:.16e}"
                           for name, _, _ in phase_points)
                + "\n")

print("Wrote velocity_profiles.dat")

# ------------------------------------------------------------
# Profile figure
# ------------------------------------------------------------

figp, (axp1, axp2) = plt.subplots(1, 2, figsize=(11, 5))

for (name, ph, q), color in zip(phase_points, ("C0", "C2", "C3")):
    axp1.plot(profiles[name], rhat, "-", lw=2, color=color,
              label=f"{name.replace('_', ' ')}  "
                    f"(t = {ph*period:.3f} s, Q = {q:.3f} L/min)")

axp1.axvline(0.0, color="0.4", lw=1)
axp1.set_xlabel("Axial velocity u_z (m/s)")
axp1.set_ylabel("r / R")
axp1.set_title("Womersley profiles imposed at the inlet")
axp1.grid(True, alpha=0.3)
axp1.legend(fontsize=8)

axp2.plot(phase * period, flow_cycle, "-", lw=2, color="C1")
axp2.axhline(0.0, color="0.4", lw=1)

for (name, ph, q), color in zip(phase_points, ("C0", "C2", "C3")):
    axp2.plot(ph * period, q, "o", ms=8, color=color,
              label=name.replace("_", " "))

axp2.set_xlabel("Time (s)")
axp2.set_ylabel("Flow rate (L/min)")
axp2.set_title("Where those profiles sit in the cycle")
axp2.grid(True, alpha=0.3)
axp2.legend(fontsize=8)

figp.tight_layout()
figp.savefig("velocity_profiles.png", dpi=150)

# ------------------------------------------------------------
# Centreline velocity and flow rate over the cycle, both rebuilt from
# the coefficients in cl.fft
#
# u_cl is what a Doppler probe on the axis would read, so this is the
# curve to compare against the digitized envelope; Q is what the
# boundary condition actually enforces.  They are not proportional --
# the profile bluntness varies through the cycle.
# ------------------------------------------------------------

expo = np.exp(2j * np.pi * np.outer(np.arange(1, NMODES + 1), phase))

u_cl_cycle = (2.0 * qmean * lpm_to_m3s / area
              + np.real((q_k * PHI[:, 0]) @ expo))

print(f"\nCentreline: min {u_cl_cycle.min():.6f}  max {u_cl_cycle.max():.6f} m/s"
      f"   (peak {u_cl_cycle.max()/u_cl_cycle.mean():.2f} x cycle mean)")

figq, axq = plt.subplots(figsize=(10, 5))

axq.plot(phase * period, u_cl_cycle, "-", lw=2, color="C3",
         label="centreline $u_z(r=0)$")

axq.axhline(0.0, color="0.4", lw=1)

for name, ph, q in phase_points:
    axq.axvline(ph * period, ls=":", lw=1, color="0.6")
    axq.annotate(name.replace("_", " "), (ph * period, axq.get_ylim()[1]),
                 fontsize=8, ha="center", va="bottom", color="0.3")

axq.set_xlabel(f"Time (s)   —   period {period:.3f} s ({freqbpm:.0f} bpm)")
axq.set_ylabel("Centreline axial velocity (m/s)")
axq.grid(True, alpha=0.3)

axq2 = axq.twinx()

axq2.plot(phase * period, flow_cycle, "--", lw=3, color="C1", alpha=0.6,
          label="flow rate $Q$")

axq2.set_ylabel("Flow rate (L/min)")

lines = axq.get_lines()[:1] + axq2.get_lines()[:1]
axq.legend(lines, [l.get_label() for l in lines], loc="upper right")

axq.set_title("Centreline velocity and flow rate from cl.fft")

figq.tight_layout()
figq.savefig("centreline_flow.png", dpi=150)


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
# Plot
# ============================================================

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 13))

for i, cyc in enumerate(cycles):
    ax1.plot(phase, cyc * qmean / cycle_mean, lw=0.8, alpha=0.5,
             label=f"beat {i+1}")

ax1.axhline(0.0, color="0.4", lw=1)
ax1.set_xlabel("Normalized phase (peak to peak)")
ax1.set_ylabel("Envelope (scaled)")
ax1.set_title("All cycles — the reverse notch drifts between phase 0.10 "
              "and 0.22,\nwhich is why averaging them removes it")
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=8, ncol=3)

ax2.plot(phase, cycle, ".", ms=4, label=f"beat {chosen+1} (digitized)")

ax2.plot(phase, cycle_reconstructed, "-", lw=2,
         label=f"{NMODES}-mode fit")

ax2.plot(phase, centreline, "--", lw=1.5,
         label="centreline rebuilt from written coefficients")

ax2.axhline(0.0, color="0.4", lw=1)
ax2.set_xlabel("Phase")
ax2.set_ylabel("Peak velocity (scaled)")
ax2.set_title(f"Measured envelope = centreline velocity — "
              f"measured {measured_bpm:.1f} bpm, "
              f"normalized to {freqbpm:.1f} bpm")
ax2.grid(True, alpha=0.3)
ax2.legend()

ax3.plot(phase, cycle, ":", lw=1.5, color="0.5",
         label="envelope shape (what the old code wrote)")

ax3.plot(phase, flow_cycle, "-", lw=2, color="C1",
         label="Q(t) after Womersley inversion")

ax3.axhline(0.0, color="0.4", lw=1)
ax3.axhline(qmean, linestyle="--", color="0.6",
            label=f"Mean = {qmean:.3f} L/min")

ax3.set_xlabel(f"Phase (end-diastole to end-diastole, "
               f"1 phase unit = {period:.3f} s)")
ax3.set_ylabel("Flow rate (L/min)")
ax3.set_title(f"Flow waveform written to cl.fft  "
              f"(peak/mean {flow_cycle.max()/qmean:.2f}, "
              f"retrograde {rev/fwd*100:.1f}%)")
ax3.grid(True, alpha=0.3)
ax3.legend()

plt.tight_layout()
plt.savefig("fft_check.png", dpi=150)
plt.show()
