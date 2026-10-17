#!/usr/bin/env python3
"""Generate tables and publication-quality figures from saved measurements."""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
import statistics

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REFERENCE = {
    "spca-small": "spca", "spca-medium": "spca", "spca-large": "spca",
    "pca-lowrank": "prcomp", "pca-flat": "prcomp", "sparse-input": "rspca-dense",
}
NAMES = {"spca": "SPCA", "rspca": "Randomized SPCA", "prcomp": "stats::prcomp",
         "rspca-dense": "Dense conversion + RSPCA", "rspca-sparse": "Sparse input + RSPCA"}
COLORS = {"spca": "#223B53", "rspca": "#007D88", "prcomp": "#223B53",
          "rspca-dense": "#223B53", "rspca-sparse": "#007D88"}


def median(rows, key):
    return statistics.median(float(row[key]) for row in rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=HERE / "results")
    args = parser.parse_args()
    with (args.results_dir / "measurements.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    groups = defaultdict(list)
    for row in rows:
        groups[(row["case"], row["method"])].append(row)
    summaries = {}
    fields = ["elapsed_seconds", "peak_rss_mib", "process_peak_rss_mib",
              "reconstruction_relative_frobenius", "fitted_energy_fraction",
              "adjusted_score_energy_fraction", "loading_zero_fraction", "iterations"]
    for (case, method), values in groups.items():
        ref = groups[(case, REFERENCE[case])]
        summary = {"case": case, "method": method, "repetitions": len(values),
                   "n": int(values[0]["n"]), "p": int(values[0]["p"]), "k": int(values[0]["k"])}
        for field in fields:
            data = [float(row[field]) for row in values]
            summary[field + "_median"] = statistics.median(data)
            summary[field + "_min"] = min(data)
            summary[field + "_max"] = max(data)
        summary["speedup_reference"] = REFERENCE[case]
        summary["speedup"] = median(ref, "elapsed_seconds") / summary["elapsed_seconds_median"]
        ref_by_rep = {row["repetition"]: float(row["elapsed_seconds"]) for row in ref}
        ratios = [ref_by_rep[row["repetition"]] / float(row["elapsed_seconds"]) for row in values]
        summary["paired_speedup_min"] = min(ratios + [summary["speedup"]])
        summary["paired_speedup_max"] = max(ratios + [summary["speedup"]])
        summary["all_converged"] = all(row["converged"] == "TRUE" for row in values)
        summary["warnings"] = sum(int(row["warnings"]) for row in values)
        # A common full-input objective uses the exact spectral scale of the
        # deterministic SPCA fit, eliminating sketch-norm scaling differences.
        if case.startswith("spca-"):
            spectral2 = median(ref, "spectral_norm_squared")
            objectives = [
                0.5 * float(row["reconstruction_error_squared"]) +
                float(row["alpha"]) * spectral2 * float(row["loading_l1"]) +
                0.5 * float(row["beta"]) * spectral2 * float(row["loading_l2_squared"])
                for row in values
            ]
            summary["full_objective_common_median"] = statistics.median(objectives)
        else:
            summary["full_objective_common_median"] = ""
        summaries[(case, method)] = summary
    with (args.results_dir / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(next(iter(summaries.values()))))
        writer.writeheader()
        writer.writerows(summaries.values())

    text = ["# Measured benchmark results", "",
            "Three isolated processes per method and workload; medians are reported. "
            "Brackets show the observed elapsed-time range. Speedup is the reference median "
            "divided by the candidate median. Full methodology is in `tools/benchmarks/README.md`.", "",
            "## Sparse loadings on approximately low-rank dense data", "",
            "| Matrix | SPCA seconds [range] | RSPCA seconds [range] | Speedup | Relative error, SPCA / RSPCA | Zero loadings |",
            "| --- | ---: | ---: | ---: | ---: | ---: |"]
    for case in ("spca-small", "spca-medium", "spca-large"):
        if (case, "spca") not in summaries:
            continue
        a, b = summaries[(case, "spca")], summaries[(case, "rspca")]
        def timing(x):
            return f'{x["elapsed_seconds_median"]:.3f} [{x["elapsed_seconds_min"]:.3f}, {x["elapsed_seconds_max"]:.3f}]'
        text.append(f'| {a["n"]:,} × {a["p"]:,} | {timing(a)} | {timing(b)} | {b["speedup"]:.2f}× | '
                    f'{a["reconstruction_relative_frobenius_median"]:.5f} / {b["reconstruction_relative_frobenius_median"]:.5f} | '
                    f'{100 * b["loading_zero_fraction_median"]:.1f}% |')
    text += ["", "All fits use rank 10, `alpha = 0.001`, `beta = 0.0001`, "
             "centering, no scaling, `tol = 1e-5`, and at most 500 iterations. "
             "RSPCA uses `o = 20`, `q = 2`. The two methods use the same relative penalties; "
             "the complete CSV also evaluates both fits using the exact deterministic spectral scale.", "",
             "## External PCA reference and a flat-spectrum limitation", "",
             "| Workload | Method | Seconds | Speedup vs prcomp | Relative error | Fitted energy / input energy |",
             "| --- | --- | ---: | ---: | ---: | ---: |"]
    for case, title in (("pca-lowrank", "Low-rank signal, 10,000 × 3,000"),
                        ("pca-flat", "Gaussian noise, 2,000 × 600")):
        for method in ("prcomp", "rspca"):
            if (case, method) not in summaries:
                continue
            value = summaries[(case, method)]
            text.append(f'| {title} | {NAMES[method]} | {value["elapsed_seconds_median"]:.3f} | '
                        f'{value["speedup"]:.2f}× | {value["reconstruction_relative_frobenius_median"]:.5f} | '
                        f'{100 * value["fitted_energy_fraction_median"]:.3f}% |')
    text += ["", "The external comparison sets `alpha = beta = 0`, so both methods target ordinary "
             "rank-10 PCA with the same preprocessing and returned scores. `prcomp(rank. = 10)` uses "
             "R's deterministic SVD; the randomized fit is approximate. This does not compare "
             "sparsity-penalized fitting with ordinary PCA.", "",
             "## Peak resident memory", "",
             "| Workload | Method | Peak RSS at fit completion (MiB) | Seconds | Relative error |",
             "| --- | --- | ---: | ---: | ---: |"]
    for case, methods, title in (("spca-large", ("spca", "rspca"), "Dense input, 10,000 × 3,000"),
                                 ("sparse-input", ("rspca-dense", "rspca-sparse"), "10% stored input, 10,000 × 3,000")):
        for method in methods:
            if (case, method) not in summaries:
                continue
            value = summaries[(case, method)]
            text.append(f'| {title} | {NAMES[method]} | {value["peak_rss_mib_median"]:.1f} | '
                        f'{value["elapsed_seconds_median"]:.3f} | {value["reconstruction_relative_frobenius_median"]:.5f} |')
    text += ["", "RSS is the whole R process high-water mark, including runtime, namespaces, input "
             "loading, a small warm-up, preprocessing, fitting, and result construction. It is read "
             "immediately after fitting and before blocked quality validation. It is neither "
             "incremental memory nor cumulative allocation. GNU time's additional whole-process "
             "peak (including validation) is also retained in the raw data.", ""]
    (args.results_dir / "tables.md").write_text("\n".join(text))

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.labelcolor": "#243744", "text.color": "#243744",
                         "xtick.color": "#243744", "ytick.color": "#243744",
                         "axes.edgecolor": "#AAB4BD", "svg.fonttype": "none"})
    png_dir = ROOT / "inst" / "figures"
    svg_dir = HERE / "figures"
    png_dir.mkdir(parents=True, exist_ok=True)
    svg_dir.mkdir(parents=True, exist_ok=True)
    def save(fig, name):
        fig.savefig(png_dir / (name + ".png"), dpi=240, bbox_inches="tight", facecolor="white")
        fig.savefig(svg_dir / (name + ".svg"), bbox_inches="tight", facecolor="white")
        plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.4), gridspec_kw={"width_ratios": [1.15, 1]})
    for method in ("spca", "rspca"):
        values = [summaries[(case, method)] for case in ("spca-small", "spca-medium", "spca-large")]
        y = np.array([v["elapsed_seconds_median"] for v in values])
        lo = y - [v["elapsed_seconds_min"] for v in values]
        hi = [v["elapsed_seconds_max"] for v in values] - y
        axes[0].errorbar(range(3), y, yerr=[lo, hi], marker="o", markersize=6,
                         linewidth=2.1, capsize=4, color=COLORS[method], label=NAMES[method])
    axes[0].set_xticks(range(3), ["2,000 × 600", "5,000 × 1,500", "10,000 × 3,000"])
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Elapsed seconds (log scale)")
    axes[0].set_xlabel("Input matrix; rank 10")
    axes[0].set_title("Sparse loadings with matched penalties", loc="left", fontweight="bold", pad=14)
    axes[0].legend(frameon=False, loc="upper left")
    axes[0].grid(axis="y", which="major", alpha=0.18)
    for index, method in enumerate(("prcomp", "rspca")):
        values = [summaries[(case, method)] for case in ("pca-lowrank", "pca-flat")]
        y = np.array([v["elapsed_seconds_median"] for v in values])
        x = np.arange(2) + (index - 0.5) * 0.33
        lo = y - [v["elapsed_seconds_min"] for v in values]
        hi = [v["elapsed_seconds_max"] for v in values] - y
        axes[1].bar(x, y, width=0.31, color=COLORS[method], label=NAMES[method],
                    yerr=[lo, hi], capsize=3)
        for xx, yy, value in zip(x, y, values):
            axes[1].text(xx, value["elapsed_seconds_max"] * 1.16,
                         f"{yy:.2f}s", ha="center", va="bottom", fontsize=9)
    axes[1].set_yscale("log")
    axes[1].set_xticks(range(2), ["Low-rank signal\n10,000 × 3,000", "Flat spectrum\n2,000 × 600"])
    axes[1].set_title("PCA against an external reference", loc="left", fontweight="bold", pad=14)
    axes[1].set_ylabel("Elapsed seconds (log scale)")
    axes[1].set_ylim(0.15, max(summaries[(case, method)]["elapsed_seconds_max"]
                             for case in ("pca-lowrank", "pca-flat")
                             for method in ("prcomp", "rspca")) * 1.85)
    axes[1].legend(frameon=False, loc="upper right")
    axes[1].grid(axis="y", which="major", alpha=0.18)
    fig.suptitle("Randomization reduces computation when a small sketch is useful", x=0.06,
                 ha="left", fontsize=14, fontweight="bold", y=1.06)
    fig.text(0.06, -0.06, "Median of 3 isolated runs; error bars show min–max. One BLAS thread. Input generation and loading excluded.", fontsize=9)
    fig.tight_layout(w_pad=3)
    save(fig, "benchmark-time")

    fig, ax = plt.subplots(figsize=(10.8, 4.3))
    cases = ("spca-small", "spca-medium", "spca-large", "pca-lowrank", "pca-flat")
    values = [summaries[(case, "rspca")] for case in cases]
    y = np.array([v["speedup"] for v in values])
    lo = y - [v["paired_speedup_min"] for v in values]
    hi = [v["paired_speedup_max"] for v in values] - y
    ax.bar(range(len(values)), y, color=["#007D88"] * 3 + ["#6A5B91"] * 2,
           width=0.62, yerr=[lo, hi], capsize=4)
    for x, (yy, value) in enumerate(zip(y, values)):
        ax.text(x, value["paired_speedup_max"] + max(y) * 0.04,
                f"{yy:.2f}×", ha="center", fontweight="bold")
    ax.axhline(1, color="#7F8A93", linewidth=1, linestyle="--")
    ax.set_xticks(range(5), ["SPCA\n2,000 × 600", "SPCA\n5,000 × 1,500", "SPCA\n10,000 × 3,000",
                            "prcomp, low rank\n10,000 × 3,000", "prcomp, flat\n2,000 × 600"])
    ax.set_ylim(0, max(v["paired_speedup_max"] for v in values) * 1.19)
    ax.set_ylabel("Reference median / RSPCA median")
    ax.set_title("Speedup depends on both the workload and the reference", loc="left", fontsize=14,
                 fontweight="bold", pad=18)
    ax.grid(axis="y", alpha=0.18)
    fig.text(0.07, -0.04, "All fits return 10 components. Error bars: paired speedup ranges over 3 runs. PCA uses zero penalties; SPCA uses matched relative penalties.", fontsize=9)
    fig.tight_layout()
    save(fig, "benchmark-speedup")

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.2), sharey=True)
    memory_panels = [("spca-large", ("spca", "rspca"), "Dense observations", ["SPCA", "RSPCA"]),
                     ("sparse-input", ("rspca-dense", "rspca-sparse"), "10% stored observations", ["Dense conversion", "Sparse operator"])]
    for ax, (case, methods, title, labels) in zip(axes, memory_panels):
        vals = [summaries[(case, method)] for method in methods]
        y = [v["peak_rss_mib_median"] for v in vals]
        ax.bar(range(2), y, width=0.55, color=["#223B53", "#007D88"])
        for x, yy in enumerate(y):
            ax.text(x, yy + 25, f"{yy:,.0f} MiB", ha="center", fontweight="bold")
        ax.set_xticks(range(2), labels)
        ax.set_title(title, loc="left", fontweight="bold", pad=14)
        ax.grid(axis="y", alpha=0.18)
        ax.margins(y=0.18)
    axes[0].set_ylabel("Peak resident memory through fit (MiB)")
    fig.suptitle("Memory use on 10,000 × 3,000 inputs", x=0.07, ha="left", y=1.04,
                 fontsize=14, fontweight="bold")
    fig.text(0.07, -0.06, "Median of 3 runs. Whole-process high-water RSS includes runtime and input storage; this is not incremental allocation.", fontsize=9)
    fig.tight_layout(w_pad=3)
    save(fig, "benchmark-memory")
    print(f"Wrote {len(summaries)} summary rows, tables, and three figures.")


if __name__ == "__main__":
    main()
