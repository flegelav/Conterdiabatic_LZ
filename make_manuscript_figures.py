"""Plot the analytic limits, field settings, and tests of the slow expansion."""
from pathlib import Path
import csv
import hashlib
import json

import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter

BASE = Path(__file__).resolve().parent
DATA = BASE / "data" / "audited_grid.csv"
OUT = BASE / "figures"
OUT.mkdir(exist_ok=True)
rows = [{key: float(value) for key, value in row.items()}
        for row in csv.DictReader(DATA.open())]
assert len(rows) == 91
assert len({(r["eps"], r["delta"]) for r in rows}) == len(rows)

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif",
    "font.size": 9,
    "axes.labelsize": 10,
    "legend.fontsize": 8,
    "axes.linewidth": .7,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "xtick.minor.size": 1.7,
    "ytick.minor.size": 1.7,
    "lines.linewidth": 1.25,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.dpi": 300,
})

mp.mp.dps = 60
K = mp.pi**2 * mp.mpf(6)**(mp.mpf(2)/3) * mp.gamma(mp.mpf(1)/3) \
    / (10 * mp.gamma(mp.mpf(2)/3)**2)

def coefficients(epsilon):
    e = mp.mpf(str(epsilon))
    a = 2 * mp.sin(mp.pi*e/6)
    b = K*e*mp.rgamma((1-e)/6)*mp.rgamma((1+e)/6)
    return float(a*a), float(2*a*b)

def save(fig, stem):
    fig.savefig(OUT / (stem + ".pdf"), bbox_inches="tight", pad_inches=.035,
                metadata={"Title": stem, "Subject": "Landau–Zener excitation and field calibration"})
    fig.savefig(OUT / (stem + ".png"), bbox_inches="tight", pad_inches=.035)
    plt.close(fig)

# Two analytic limits and the calculated field settings between them.
from numerics.independent_check import endpoint_quadrature, inverse_settings

with (BASE / "data" / "calibration_tolerance.csv").open() as stream:
    settings = list(csv.DictReader(stream))
delta = np.array([float(row["delta"]) for row in settings])
error_percent = np.array([100 * float(row["propagated_epsilon"]) for row in settings])
target = float(settings[0]["target_ratio"])
if any(float(row["target_ratio"]) != target for row in settings):
    raise ValueError("Field settings must correspond to a single target")
leading, shift = inverse_settings(target)
fast = 2 * np.arcsin(np.sqrt(target)) / np.pi
curve_rates = np.geomspace(min(delta), max(delta), 257)
curve = 100 * np.sqrt(target / np.array([endpoint_quadrature(d) for d in curve_rates]))
small_rates = np.geomspace(min(delta), .1, 200)
small_prediction = 100 * (leading - shift * small_rates**(2/3))
eta = np.linspace(0, 2, 501)
slow_color, fast_color, ink = "#146884", "#A44824", "#26292E"

fig, (response, tolerance) = plt.subplots(1, 2, figsize=(7.05, 2.94),
                                         gridspec_kw={"width_ratios": [1, 1.51]})
fig.subplots_adjust(left=.091, right=.989, bottom=.188, top=.94, wspace=.36)
response.plot(eta, 4*np.sin(np.pi*(1-eta)/6)**2, color=slow_color, lw=1.9,
              label="Slow sweep")
response.plot(eta, np.sin(np.pi*(1-eta)/2)**2, "--", color=fast_color, lw=1.7,
              label="Fast sweep")
response.set(xlim=(0, 2), ylim=(0, 1.04), xticks=[0, .5, 1, 1.5, 2],
             yticks=[0,.25,.5,.75,1], ylabel=r"Surviving fraction $P_\eta/P_{\rm LZ}$",
             xlabel=r"Gain $\eta$ (applied / ideal amplitude)")
response.legend(loc="upper center", frameon=False, labelspacing=.2, handlelength=1.6)
response.text(.04,.08,"(a)",transform=response.transAxes)
response.grid(axis="y",color=".9",lw=.5,zorder=0)

tolerance.set_xscale("log")
tolerance.axhline(100*leading, ls=":", color=slow_color, lw=1.4)
tolerance.axhline(100*fast, ls="--", color=fast_color, lw=1.25)
tolerance.plot(curve_rates, curve, color="#4F5685", lw=1.6, label="From curvature")
tolerance.plot(delta, error_percent, "o", color=ink, mfc="white", ms=3.15, mew=.75,
               label="Propagation")
tolerance.set(xlim=(.00075, 1.5e5), ylim=(6.13, 9.91),
              yticks=[6.5,7.5,8.5,9.5], xticks=[.001,.1,10,1000,100000],
              xlabel=r"Sweep rate $\delta$", ylabel="Gain tolerance (%)")
tolerance.text(.976,.952,f"Slow  {100*leading:.2f}%",ha="right",va="center",
               transform=tolerance.transAxes,color=slow_color,fontsize=8)
tolerance.text(.035,.094,f"Fast  {100*fast:.2f}%",transform=tolerance.transAxes,
               color=fast_color,fontsize=8)
tolerance.text(.025,.80,"(b)",transform=tolerance.transAxes)
tolerance.legend(loc="upper right", bbox_to_anchor=(1,.872), frameon=False,
                 fontsize=7.6, handlelength=1.8, borderaxespad=.2, labelspacing=.22)
mark = np.flatnonzero(delta == 2)[0]
tolerance.annotate(f"{error_percent[mark]:.2f}%", xy=(delta[mark],error_percent[mark]),
                   xytext=(12,8.47),fontsize=8.4,color=ink,
                   arrowprops=dict(arrowstyle="-",lw=.66,color=ink,shrinkA=0,shrinkB=3))
zoom = tolerance.inset_axes([.115,.24,.295,.27])
zoom.set_xscale("log")
zoom.plot(small_rates, small_prediction, lw=1.25, color=slow_color)
zoom.axhline(100*leading,color=slow_color,ls=":",lw=.75)
mask = delta <= .1
zoom.plot(delta[mask], error_percent[mask], "o", color=ink, mfc="white", ms=2.5,mew=.65)
zoom.set(xlim=(.00075,.14),ylim=(9.22,9.585),xticks=[.001,.1],
         xticklabels=[r"$10^{-3}$",r"$10^{-1}$"],yticks=[9.3,9.5])
zoom.tick_params(labelsize=6.8,pad=1.2,length=2)
zoom.set_title(r"Small-$\delta$ expansion", fontsize=7.2, pad=2)
for spine in zoom.spines.values():
    spine.set_linewidth(.55)
save(fig,"calibration_ratio")

# Fixed-epsilon comparison on the original grid.
half = sorted([r for r in rows if r["eps"] == .5], key=lambda r: r["delta"])
x = np.array([r["delta"]**(2/3) for r in half])
C = np.array([r["C"] for r in half])
C0, dC = coefficients(.5)
grid = np.linspace(0, .5**(2/3), 250)
fig, ax = plt.subplots(figsize=(3.45, 2.65))
ax.plot(grid, C0 + dC*grid, color="#0072B2", label="Two-term law", zorder=2)
ax.axhline(C0, color="#555555", ls="--", label="Leading law", zorder=1)
ax.axhline(.5*np.exp(-1/3), color="#B4372D", ls=":", lw=1.5,
           label="Cardoso (2024)", zorder=1)
ax.plot(x, C, "o", color="black", ms=3.0, markeredgewidth=.5,
        label="Propagation", zorder=3)
handles, labels = ax.get_legend_handles_labels()
ax.legend([handles[i] for i in [3, 0, 1, 2]], [labels[i] for i in [3, 0, 1, 2]],
          loc="center left", bbox_to_anchor=(.025, .67), frameon=False,
          borderaxespad=0, handlelength=2.1, labelspacing=.30)
ax.text(.97, .96, r"$\varepsilon=0.5$", transform=ax.transAxes,
        ha="right", va="top", fontsize=9)
ax.set(xlabel=r"$\delta^{2/3}$", ylabel=r"Surviving fraction $P_\eta/P_{\rm LZ}$",
       xlim=(-.012, .66), ylim=(.262, .369))
ax.set_xticks([0, .2, .4, .6])
ax.set_yticks([.27, .30, .33, .36])
fig.tight_layout(pad=.35)
save(fig, "prefactor_comparison")

# Supplemental figure: relative probability errors across representative epsilon values.
# No slope fitting or numerical discrepancy floor is shown: neither is needed
# to display the improvement, and the proof only asserts an O(delta) remainder.
selected = [(.001, "#0072B2", "o"), (.5, "#D55E00", "s"),
            (.95, "#009E73", "^")]
fig, axes = plt.subplots(1, 2, figsize=(6.95, 2.65), sharex=True, sharey=True)
for epsilon, color, marker in selected:
    rr = sorted([r for r in rows if r["eps"] == epsilon], key=lambda r: r["delta"])
    ds = [r["delta"] for r in rr]
    for ax, key in zip(axes, ["error_leading", "error_two_term"]):
        ax.loglog(ds, [r[key] for r in rr], color=color, marker=marker,
                  markersize=3.2, markerfacecolor="white", markeredgewidth=.8,
                  label=fr"$\varepsilon={epsilon:g}$")
for i, (ax, title) in enumerate(zip(axes, ["Leading law", "Two-term law"])):
    ax.text(.04, .95, f"({chr(97+i)}) {title}", transform=ax.transAxes,
            ha="left", va="top", fontsize=9)
    ax.set(xlim=(.0008, .65), ylim=(4e-7, .25), xlabel=r"$\delta$")
    ax.grid(which="major", color=".91", lw=.5)
    ax.xaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_minor_formatter(NullFormatter())
axes[0].set_ylabel("Relative probability error")
axes[1].legend(loc="lower right", frameon=False, labelspacing=.35,
               handlelength=2, borderaxespad=.55)
fig.tight_layout(pad=.35, w_pad=.8)
save(fig, "relative_errors")

summary = {
    "input_sha256": hashlib.sha256(DATA.read_bytes()).hexdigest(),
    "grid_points": len(rows),
    "grid_runs_in_source_audit": 182,
    "epsilon_values": sorted({r["eps"] for r in rows}),
    "delta_values": sorted({r["delta"] for r in rows}),
    "plotted_epsilon_values_supplemental_figure": [x[0] for x in selected],
    "epsilon_half_table": [{k: r[k] for k in
        ["delta", "C", "error_leading", "error_two_term", "error_cardoso"]}
        for r in reversed(half) if r["delta"] in [.1, .01, .001]],
    "target_for_field_settings": target,
    "propagated_field_settings": len(settings),
    "field_settings_sha256": hashlib.sha256((BASE / "data/calibration_tolerance.csv").read_bytes()).hexdigest(),
    "curvature_quadrature_points_in_fig1b": len(curve_rates),
    "curvature_check_sha256": hashlib.sha256((BASE / "data/curvature_tolerance.csv").read_bytes()).hexdigest(),
    "scope": "The curve in Fig. 1(b) uses curvature quadrature at logarithmically spaced rates. The circles give propagated thresholds, horizontal lines give limiting thresholds, and the inset shows the slow rate correction."
}
(BASE / "figure_data_summary.json").write_text(json.dumps(summary, indent=2)+"\n")
print(json.dumps(summary, indent=2))
