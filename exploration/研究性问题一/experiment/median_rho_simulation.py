"""Monte Carlo comparison of MedRank and exact Spearman-rho aggregation.

profile[i, a] is the 1-based rank of alternative a assigned by voter i.
Each row is a permutation. Following the uniform-input experiment in Section 5
of Bednay, Fleiner and Tasnadi (2026), the default experiment now considers
10 through 20 alternatives, 3 through 50 voters, and 2500 profiles per point.

MedRank uses upper medians and fixed alternative-index tie breaking. The exact
rho optimum is obtained by Borda rank-sum sorting (Proposition 7.1), which avoids
a minimum-cost assignment solver. Batch computation changes neither algorithm.
Sample maxima certify lower bounds, not upper bounds, for the worst-case ratio.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np


def ranks_from_keys(keys: np.ndarray) -> np.ndarray:
    """Convert the last axis of ascending keys to ranks, breaking ties by index."""
    order = np.argsort(keys, axis=-1, kind="stable")
    result = np.empty_like(order)
    positions = np.broadcast_to(np.arange(1, order.shape[-1] + 1), order.shape)
    np.put_along_axis(result, order, positions, axis=-1)
    return result


def random_profile(rng: np.random.Generator, voters: int, alternatives: int) -> np.ndarray:
    """Sample voters independently and uniformly from all m! rank vectors."""
    return np.argsort(rng.random((voters, alternatives)), axis=1) + 1


def medrank(profile: np.ndarray) -> np.ndarray:
    """Rank by upper median; break equal median ranks by alternative index.

    For an even number of voters, the upper median is the smallest rank at
    which STRICTLY more than half the voters rank the alternative at or above
    that position.  This is the MedRank convention used in the cited paper.
    """
    voters = profile.shape[0]
    median_ranks = np.sort(profile, axis=0)[voters // 2]
    return ranks_from_keys(median_ranks)


def rho_cost(profile: np.ndarray, ranking: np.ndarray) -> int:
    """Unnormalized total squared rank difference over all voters and objects."""
    delta = profile.astype(np.int64) - ranking.astype(np.int64)
    return int(np.sum(delta * delta))


def rho_optimum(profile: np.ndarray) -> tuple[np.ndarray, int]:
    """Return an exact rho optimum by sorting rank sums (Borda).

    Expanding the square leaves a constant minus twice the scalar product of
    the output ranks and the input rank sums. The rearrangement inequality
    therefore makes their common ascending order optimal. Equal sums are
    resolved by alternative index; every such ordering has the same cost.
    """
    optimum = ranks_from_keys(profile.sum(axis=0, dtype=np.int64))
    return optimum, rho_cost(profile, optimum)


def evaluate_batch(profiles: np.ndarray) -> tuple[np.ndarray, ...]:
    """Evaluate (samples, voters, alternatives) profiles without Python loops."""
    voters = profiles.shape[1]
    medians = np.partition(profiles, voters // 2, axis=1)[:, voters // 2, :]
    med = ranks_from_keys(medians)
    optimal = ranks_from_keys(profiles.sum(axis=1, dtype=np.int64))
    med_delta = profiles - med[:, None, :]
    opt_delta = profiles - optimal[:, None, :]
    med_costs = np.sum(med_delta * med_delta, axis=(1, 2), dtype=np.int64)
    opt_costs = np.sum(opt_delta * opt_delta, axis=(1, 2), dtype=np.int64)
    return med, optimal, med_costs, opt_costs


def standard_error(values: np.ndarray) -> float:
    """Sample standard error; one observation cannot estimate variance."""
    if values.size < 2:
        return float("nan")
    return float(np.std(values, ddof=1) / np.sqrt(values.size))


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_example(path: Path, example: dict | None) -> None:
    """Always replace the witness file, including when all optima are zero."""
    if example is None:
        path.write_text("No profile with a positive optimum was sampled.\n", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8") as file:
        file.write(f"alternatives: {example['alternatives']}\n"
                   f"voters: {example['voters']}\n"
                   f"base seed: {example['seed']}\n"
                   f"seed entropy: [{example['seed']}, {example['alternatives']}]\n"
                   f"ratio: {example['ratio']:.12f}\n"
                   f"MedRank cost: {example['medrank_cost']}\n"
                   f"optimum cost: {example['optimum_cost']}\n"
                   "Each row is a voter; columns are alternatives 1, ..., m.\nprofile:\n")
        np.savetxt(file, example["profile"], fmt="%d")
        file.write(f"MedRank rank vector: {example['medrank']}\n"
                   f"optimal rank vector: {example['optimal']}\n")


def run_experiment(args: argparse.Namespace, alternatives: int) -> tuple[list[dict], dict, dict | None]:
    # The stream for each object count is independent of the requested list order.
    rng = np.random.default_rng(np.random.SeedSequence([args.seed, alternatives]))
    maximum_distance = (alternatives**3 - alternatives) // 3
    rows: list[dict] = []
    all_ratios: list[np.ndarray] = []
    largest: dict | None = None
    for voters in range(args.min_voters, args.max_voters + 1):
        med_costs = np.empty(args.samples, dtype=np.int64)
        opt_costs = np.empty(args.samples, dtype=np.int64)
        for start in range(0, args.samples, args.batch_size):
            count = min(args.batch_size, args.samples - start)
            profiles = np.argsort(rng.random((count, voters, alternatives)), axis=2) + 1
            med, optimal, med_batch, opt_batch = evaluate_batch(profiles)
            assert np.all(med_batch >= opt_batch)
            assert np.all(med_batch[opt_batch == 0] == 0)
            med_costs[start:start + count] = med_batch
            opt_costs[start:start + count] = opt_batch
            ratios = np.full(count, np.nan)
            np.divide(med_batch, opt_batch, out=ratios, where=opt_batch > 0)
            if np.any(opt_batch > 0):
                index = int(np.nanargmax(ratios))
                if largest is None or ratios[index] > largest["ratio"]:
                    largest = {
                        "alternatives": alternatives, "voters": voters, "seed": args.seed,
                        "ratio": float(ratios[index]),
                        "medrank_cost": int(med_batch[index]), "optimum_cost": int(opt_batch[index]),
                        "profile": profiles[index].tolist(),
                        "medrank": med[index].tolist(), "optimal": optimal[index].tolist(),
                    }
        positive = opt_costs > 0
        ratios = med_costs[positive] / opt_costs[positive]
        all_ratios.append(ratios)
        mean_ratio = float(np.mean(ratios)) if ratios.size else float("nan")
        ratio_se = standard_error(ratios)
        quantiles = np.quantile(ratios, [0.5, 0.9, 0.95, 0.99]) if ratios.size else [float("nan")] * 4
        normalizer = voters * maximum_distance
        gaps = (med_costs - opt_costs) / normalizer
        gap_mean = float(np.mean(gaps))
        gap_se = standard_error(gaps)
        row = {
            "alternatives": alternatives, "voters": voters, "samples": args.samples,
            "zero_optimum_profiles": int(np.count_nonzero(~positive)),
            "mean_medrank_normalized": float(np.mean(med_costs)) / normalizer,
            "mean_optimum_normalized": float(np.mean(opt_costs)) / normalizer,
            "mean_profile_ratio": mean_ratio,
            "mean_profile_ratio_se": ratio_se,
            "mean_profile_ratio_ci95_low": mean_ratio - 1.96 * ratio_se,
            "mean_profile_ratio_ci95_high": mean_ratio + 1.96 * ratio_se,
            "ratio_p50": float(quantiles[0]), "ratio_p90": float(quantiles[1]),
            "ratio_p95": float(quantiles[2]), "ratio_p99": float(quantiles[3]),
            "max_sampled_ratio": float(np.max(ratios)) if ratios.size else float("nan"),
            "mean_normalized_gap": gap_mean, "mean_normalized_gap_se": gap_se,
            "mean_normalized_gap_ci95_low": gap_mean - 1.96 * gap_se,
            "mean_normalized_gap_ci95_high": gap_mean + 1.96 * gap_se,
        }
        rows.append(row)
        if voters == args.min_voters or voters == args.max_voters or voters % 10 == 0:
            print(f"alternatives={alternatives:2d} voters={voters:3d} "
                  f"mean ratio={mean_ratio:.6f} max ratio={row['max_sampled_ratio']:.6f}", flush=True)
    joined_ratios = np.concatenate(all_ratios)
    finite_mean_ratios = [row["mean_profile_ratio"] for row in rows
                          if np.isfinite(row["mean_profile_ratio"])]
    comparison = {
        "alternatives": alternatives, "total_profiles": args.samples * len(rows),
        "mean_profile_ratio_all_voters": float(np.mean(joined_ratios)) if joined_ratios.size else float("nan"),
        "min_mean_profile_ratio": min(finite_mean_ratios, default=float("nan")),
        "max_mean_profile_ratio": max(finite_mean_ratios, default=float("nan")),
        "overall_ratio_p95": float(np.quantile(joined_ratios, 0.95)) if joined_ratios.size else float("nan"),
        "overall_ratio_p99": float(np.quantile(joined_ratios, 0.99)) if joined_ratios.size else float("nan"),
        "max_sampled_ratio": largest["ratio"] if largest else float("nan"),
        "largest_example_voters": largest["voters"] if largest else "",
    }
    folder = args.out / f"alternatives_{alternatives}"
    folder.mkdir(parents=True, exist_ok=True)
    write_csv(folder / "summary.csv", rows)
    write_example(folder / "largest_sampled_ratio.txt", largest)
    (folder / "largest_sampled_ratio.json").write_text(json.dumps(largest, indent=2), encoding="utf-8")
    return rows, comparison, largest


def save_plots(out: Path, groups: dict[int, list[dict]], comparisons: list[dict], samples: int) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for alternatives, rows in groups.items():
        folder = out / f"alternatives_{alternatives}"
        x = [row["voters"] for row in rows]
        fig, ax = plt.subplots(figsize=(8, 4.8))
        ax.plot(x, [row["mean_medrank_normalized"] for row in rows], label="MedRank")
        ax.plot(x, [row["mean_optimum_normalized"] for row in rows], label="Exact rho optimum (Borda)")
        ax.set(xlabel="Number of voters", ylabel="Mean normalized rho distance",
               title=f"{alternatives} alternatives; {samples} profiles per point")
        ax.grid(alpha=0.25)
        ax.legend()
        fig.tight_layout()
        fig.savefig(folder / "mean_distances.png", dpi=180)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8, 4.8))
        ax.plot(x, [row["mean_profile_ratio"] for row in rows], label="Mean profile ratio")
        ax.fill_between(x, [row["mean_profile_ratio_ci95_low"] for row in rows],
                        [row["mean_profile_ratio_ci95_high"] for row in rows], alpha=0.2,
                        label="Approximate 95% interval for mean")
        ax.plot(x, [row["max_sampled_ratio"] for row in rows], label="Largest sampled ratio")
        ax.set(xlabel="Number of voters", ylabel="MedRank cost / optimum cost",
               title=f"Sampled ratios: {alternatives} alternatives")
        ax.grid(alpha=0.25)
        ax.legend()
        fig.tight_layout()
        fig.savefig(folder / "ratios.png", dpi=180)
        plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    counts = list(groups)
    available_voters = {row["voters"] for row in next(iter(groups.values()))}
    selected_voters = [k for k in [3, 4, 10, 20, 50] if k in available_voters]
    if not selected_voters:
        selected_voters = [min(available_voters)]
    for voters in selected_voters:
        points = [next(row for row in groups[m] if row["voters"] == voters) for m in counts]
        axes[0].errorbar(counts, [row["mean_profile_ratio"] for row in points],
                         yerr=[1.96 * row["mean_profile_ratio_se"] for row in points],
                         marker="o", capsize=2, label=f"{voters} voters")
    axes[0].set(xlabel="Number of alternatives", ylabel="Mean profile ratio",
                title="Mean ratios with approximate 95% intervals")
    axes[0].legend(fontsize=8)
    axes[1].plot(counts, [row["max_sampled_ratio"] for row in comparisons],
                 marker="o", label="Largest ratio over all voter counts")
    axes[1].plot(counts, [row["overall_ratio_p99"] for row in comparisons],
                 marker="o", label="99th percentile over all voter counts")
    axes[1].set(xlabel="Number of alternatives", ylabel="Sampled profile ratio",
                title="Observed upper tail of random inputs")
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=0.25)
        ax.set_xticks(counts)
    fig.tight_layout()
    fig.savefig(out / "comparison_by_alternatives.png", dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alternatives", type=int, nargs="+", default=list(range(10, 21)),
                        help="object counts to simulate; default: every integer from 10 to 20")
    parser.add_argument("--min-voters", type=int, default=3)
    parser.add_argument("--max-voters", type=int, default=50)
    parser.add_argument("--samples", type=int, default=2500,
                        help="independent preference profiles per voter count")
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--batch-size", type=int, default=250)
    parser.add_argument("--out", type=Path,
                        default=Path(__file__).resolve().parent / "median_rho_results_10_20")
    args = parser.parse_args()
    if (min(args.alternatives) < 2 or args.min_voters < 1 or args.max_voters < args.min_voters
            or args.samples < 1 or args.batch_size < 1 or args.seed < 0):
        parser.error("Require alternatives >= 2, voters >= 1, max >= min, samples/batch >= 1, seed >= 0.")
    args.alternatives = sorted(set(args.alternatives))

    args.out.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    config = {
        "alternatives": args.alternatives, "min_voters": args.min_voters,
        "max_voters": args.max_voters, "samples_per_point": args.samples,
        "seed": args.seed, "seed_entropy_per_alternative_count": "[seed, alternatives]",
        "numpy_bit_generator": "PCG64", "batch_size": args.batch_size,
        "sampling": "independent uniform random permutations",
        "medrank": "upper median, then alternative-index tie breaking",
        "optimum": "exact rho minimum via stable Borda rank-sum sorting",
        "normalizer": "voters * (alternatives**3 - alternatives) / 3",
        "ci95": "approximate pointwise mean +/- 1.96 * sample standard error; not simultaneous",
        "python_version": platform.python_version(), "numpy_version": np.__version__,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    groups: dict[int, list[dict]] = {}
    comparisons: list[dict] = []
    largest: dict | None = None
    for alternatives in args.alternatives:
        rows, comparison, example = run_experiment(args, alternatives)
        groups[alternatives] = rows
        comparisons.append(comparison)
        if example is not None and (largest is None or example["ratio"] > largest["ratio"]):
            largest = example
    write_csv(args.out / "summary.csv", [row for rows in groups.values() for row in rows])
    write_csv(args.out / "comparison.csv", comparisons)
    write_example(args.out / "largest_sampled_ratio.txt", largest)
    (args.out / "largest_sampled_ratio.json").write_text(json.dumps(largest, indent=2), encoding="utf-8")
    save_plots(args.out, groups, comparisons, args.samples)
    import matplotlib
    config.update({"matplotlib_version": matplotlib.__version__,
                   "total_profiles": args.samples * (args.max_voters - args.min_voters + 1) * len(groups),
                   "elapsed_seconds": time.perf_counter() - started,
                   "finished_utc": datetime.now(timezone.utc).isoformat()})
    (args.out / "run_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
    print(f"Completed {config['total_profiles']:,} profiles in {config['elapsed_seconds']:.1f}s; "
          f"results: {args.out}", flush=True)


if __name__ == "__main__":
    main()
