#!/usr/bin/env python3
"""Python audit for Micro-CT spicule-arrangement morphometric analyses.

This script mirrors the public R workflow as closely as possible using the
bundled Python runtime available in Codex. It is intended as a reproducible
audit when R/vegan is unavailable. Final manuscript p-values can still be
regenerated with the R script and vegan before submission.
"""

from __future__ import annotations

import argparse
import itertools
import math
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_SEED = 20260707
DEFAULT_PERMUTATIONS = 9999
RUN_LABEL = "morphometric_analysis_20260806_R002_R003_D"


def normal_two_sided_p(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2.0))


def chi_square_sf_df3(x: float) -> float:
    # Survival function for chi-square distribution with df = 3.
    # CDF = erf(sqrt(x/2)) - sqrt(2x/pi) exp(-x/2).
    if x <= 0:
        return 1.0
    return math.erfc(math.sqrt(x / 2.0)) + math.sqrt(2.0 * x / math.pi) * math.exp(-x / 2.0)


def holm_adjust(p_values: np.ndarray) -> np.ndarray:
    p_values = np.asarray(p_values, dtype=float)
    order = np.argsort(p_values)
    adjusted = np.empty_like(p_values)
    running_max = 0.0
    m = len(p_values)
    for rank, idx in enumerate(order):
        value = min(1.0, (m - rank) * p_values[idx])
        running_max = max(running_max, value)
        adjusted[idx] = running_max
    return adjusted


def average_ranks(values: np.ndarray) -> tuple[np.ndarray, float]:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    tie_sum = 0.0
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and values[order[end]] == values[order[start]]:
            end += 1
        avg_rank = (start + 1 + end) / 2.0
        ranks[order[start:end]] = avg_rank
        tie_len = end - start
        if tie_len > 1:
            tie_sum += tie_len**3 - tie_len
        start = end
    return ranks, tie_sum


def standardize(df: pd.DataFrame) -> pd.DataFrame:
    means = df.mean(axis=0)
    sds = df.std(axis=0, ddof=1)
    if (sds == 0).any():
        zero_cols = ", ".join(sds.index[sds == 0])
        raise ValueError(f"Cannot standardize zero-variance columns: {zero_cols}")
    return (df - means) / sds


def pca_from_scaled(x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    u, s, vt = np.linalg.svd(x, full_matrices=False)
    scores = u * s
    loadings = vt.T
    sdev = s / math.sqrt(x.shape[0] - 1)
    variance = (sdev**2) / np.sum(sdev**2)
    return scores, loadings, sdev, variance


def permanova_stats(x: np.ndarray, labels: np.ndarray) -> dict[str, float]:
    unique = np.unique(labels)
    n = x.shape[0]
    k = len(unique)
    overall = x.mean(axis=0)
    total_ss = float(np.sum((x - overall) ** 2))
    between_ss = 0.0
    for group in unique:
        xg = x[labels == group]
        between_ss += xg.shape[0] * float(np.sum((xg.mean(axis=0) - overall) ** 2))
    within_ss = total_ss - between_ss
    df_between = k - 1
    df_within = n - k
    f_value = (between_ss / df_between) / (within_ss / df_within)
    return {
        "Df": float(df_between),
        "Df_residual": float(df_within),
        "SumOfSqs": between_ss,
        "SumOfSqs_residual": within_ss,
        "R2": between_ss / total_ss,
        "F": f_value,
    }


def permanova_test(
    x: np.ndarray,
    labels: np.ndarray,
    permutations: int,
    rng: np.random.Generator,
) -> dict[str, float]:
    observed = permanova_stats(x, labels)
    f_observed = observed["F"]
    exceed = 0
    labels = np.asarray(labels)
    for _ in range(permutations):
        permuted = rng.permutation(labels)
        if permanova_stats(x, permuted)["F"] >= f_observed - 1e-12:
            exceed += 1
    observed["P"] = (exceed + 1) / (permutations + 1)
    return observed


def geometric_median(x: np.ndarray, tolerance: float = 1e-9, max_iter: int = 1000) -> np.ndarray:
    y = np.median(x, axis=0)
    for _ in range(max_iter):
        distances = np.linalg.norm(x - y, axis=1)
        if np.any(distances < tolerance):
            return x[np.argmin(distances)]
        weights = 1.0 / distances
        y_next = np.sum(x * weights[:, None], axis=0) / np.sum(weights)
        if np.linalg.norm(y_next - y) < tolerance:
            return y_next
        y = y_next
    return y


def permdisp_distances(x: np.ndarray, labels: np.ndarray) -> np.ndarray:
    distances = np.empty(x.shape[0], dtype=float)
    for group in np.unique(labels):
        mask = labels == group
        median = geometric_median(x[mask])
        distances[mask] = np.linalg.norm(x[mask] - median, axis=1)
    return distances


def one_way_anova_f(values: np.ndarray, labels: np.ndarray) -> dict[str, float]:
    unique = np.unique(labels)
    n = len(values)
    k = len(unique)
    overall = float(np.mean(values))
    between = 0.0
    within = 0.0
    for group in unique:
        vg = values[labels == group]
        mean_g = float(np.mean(vg))
        between += len(vg) * (mean_g - overall) ** 2
        within += float(np.sum((vg - mean_g) ** 2))
    df_between = k - 1
    df_within = n - k
    f_value = (between / df_between) / (within / df_within)
    return {
        "Df": float(df_between),
        "Df_residual": float(df_within),
        "SumOfSqs": between,
        "SumOfSqs_residual": within,
        "F": f_value,
    }


def permdisp_test(
    x: np.ndarray,
    labels: np.ndarray,
    permutations: int,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, dict[str, float]]:
    observed_distances = permdisp_distances(x, labels)
    observed = one_way_anova_f(observed_distances, labels)
    exceed = 0
    for _ in range(permutations):
        permuted = rng.permutation(labels)
        perm_distances = permdisp_distances(x, permuted)
        if one_way_anova_f(perm_distances, permuted)["F"] >= observed["F"] - 1e-12:
            exceed += 1
    observed["P"] = (exceed + 1) / (permutations + 1)
    distance_df = pd.DataFrame({"Network": labels, "DistanceToMedian": observed_distances})
    return distance_df, observed


def kruskal_wallis(values: np.ndarray, labels: np.ndarray) -> dict[str, float]:
    ranks, tie_sum = average_ranks(values)
    n = len(values)
    unique = np.unique(labels)
    rank_term = 0.0
    for group in unique:
        rg = ranks[labels == group]
        rank_term += float(np.sum(rg) ** 2) / len(rg)
    h = (12.0 / (n * (n + 1.0))) * rank_term - 3.0 * (n + 1.0)
    tie_correction = 1.0 - tie_sum / (n**3 - n)
    h_corrected = h / tie_correction
    return {
        "Statistic": h_corrected,
        "Df": float(len(unique) - 1),
        "P": chi_square_sf_df3(h_corrected),
    }


def dunn_test(values: np.ndarray, labels: np.ndarray, variable: str) -> pd.DataFrame:
    ranks, tie_sum = average_ranks(values)
    n = len(values)
    unique = np.unique(labels)
    variance = n * (n + 1.0) / 12.0 - tie_sum / (12.0 * (n - 1.0))
    rows = []
    for group1, group2 in itertools.combinations(unique, 2):
        mask1 = labels == group1
        mask2 = labels == group2
        z = (np.mean(ranks[mask1]) - np.mean(ranks[mask2])) / math.sqrt(
            variance * (1.0 / np.sum(mask1) + 1.0 / np.sum(mask2))
        )
        rows.append(
            {
                "Variable": variable,
                "Comparison": f"{group1} - {group2}",
                "Z": z,
                "P": normal_two_sided_p(z),
            }
        )
    out = pd.DataFrame(rows)
    out["P_adjusted_holm"] = holm_adjust(out["P"].to_numpy())
    return out


def format_table(df: pd.DataFrame, digits: int = 6) -> str:
    return df.to_string(index=False, float_format=lambda x: f"{x:.{digits}g}")


def run(input_csv: Path, output_dir: Path, permutations: int, seed: int) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(input_csv)
    required = {"Sample", "Network"}
    missing = required.difference(raw.columns)
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(sorted(missing))}")

    response_columns = [c for c in raw.columns if c not in required]
    responses = raw[response_columns].apply(pd.to_numeric, errors="coerce")
    keep = responses.notna().all(axis=1) & raw["Network"].notna()
    raw = raw.loc[keep].copy()
    responses = responses.loc[keep].copy()
    labels = raw["Network"].astype(str).to_numpy()
    samples = raw["Sample"].astype(str).to_numpy()
    scaled = standardize(responses)
    x = scaled.to_numpy(dtype=float)

    rng = np.random.default_rng(seed)

    scores, loadings, sdev, variance = pca_from_scaled(x)
    pc_names = [f"PC{i}" for i in range(1, len(variance) + 1)]
    pca_scores = pd.DataFrame(scores, columns=pc_names)
    pca_scores.insert(0, "Network", labels)
    pca_scores.insert(0, "Sample", samples)
    pca_loadings = pd.DataFrame(loadings, columns=pc_names)
    pca_loadings.insert(0, "Variable", response_columns)
    pca_variance = pd.DataFrame(
        {
            "PC": pc_names,
            "StandardDeviation": sdev,
            "ProportionVariance": variance,
            "CumulativeVariance": np.cumsum(variance),
        }
    )

    global_permanova = pd.DataFrame([permanova_test(x, labels, permutations, rng)])

    pairwise_rows = []
    for group1, group2 in itertools.combinations(np.unique(labels), 2):
        mask = (labels == group1) | (labels == group2)
        stats = permanova_test(x[mask], labels[mask], permutations, rng)
        stats.update({"Group1": group1, "Group2": group2})
        pairwise_rows.append(stats)
    pairwise = pd.DataFrame(pairwise_rows)
    pairwise["P_adjusted_holm"] = holm_adjust(pairwise["P"].to_numpy())
    pairwise = pairwise[
        ["Group1", "Group2", "Df", "Df_residual", "SumOfSqs", "SumOfSqs_residual", "R2", "F", "P", "P_adjusted_holm"]
    ]

    permdisp_dist, permdisp_stats = permdisp_test(x, labels, permutations, rng)
    permdisp_dist.insert(0, "Sample", samples)
    permdisp_global = pd.DataFrame([permdisp_stats])

    kruskal_rows = []
    dunn_rows = []
    for variable in response_columns:
        values = responses[variable].to_numpy(dtype=float)
        kw = kruskal_wallis(values, labels)
        kw["Variable"] = variable
        kruskal_rows.append(kw)
        dunn_rows.append(dunn_test(values, labels, variable))
    kruskal = pd.DataFrame(kruskal_rows)[["Variable", "Statistic", "Df", "P"]]
    kruskal["P_adjusted_holm"] = holm_adjust(kruskal["P"].to_numpy())
    dunn = pd.concat(dunn_rows, ignore_index=True)

    pca_scores.to_csv(output_dir / "pca_scores.csv", index=False)
    pca_loadings.to_csv(output_dir / "pca_loadings.csv", index=False)
    pca_variance.to_csv(output_dir / "pca_variance.csv", index=False)
    global_permanova.to_csv(output_dir / "permanova_global.csv", index=False)
    pairwise.to_csv(output_dir / "pairwise_permanova.csv", index=False)
    permdisp_dist.to_csv(output_dir / "permdisp_distances_to_median.csv", index=False)
    permdisp_global.to_csv(output_dir / "permdisp_global.csv", index=False)
    kruskal.to_csv(output_dir / "kruskal_wallis_by_variable.csv", index=False)
    dunn.to_csv(output_dir / "dunn_tests_by_variable.csv", index=False)

    group_sizes = pd.Series(labels).value_counts().sort_index()
    with (output_dir / "summary.txt").open("w", encoding="utf-8") as handle:
        handle.write(f"Input file: {input_csv}\n")
        handle.write(f"Output directory: {output_dir}\n")
        handle.write(f"Seed: {seed}\n")
        handle.write(f"Permutations: {permutations}\n")
        handle.write(f"Samples included: {len(samples)}\n")
        handle.write(f"Variables included: {len(response_columns)}\n")
        handle.write("\nNetwork group sizes:\n")
        handle.write(group_sizes.to_string())
        handle.write("\n\nPCA variance summary:\n")
        handle.write(format_table(pca_variance.head(12)))
        handle.write("\n\nPERMANOVA global:\n")
        handle.write(format_table(global_permanova))
        handle.write("\n\nPairwise PERMANOVA:\n")
        handle.write(format_table(pairwise))
        handle.write("\n\nPERMDISP global:\n")
        handle.write(format_table(permdisp_global))
        handle.write("\n\nKruskal-Wallis:\n")
        handle.write(format_table(kruskal))
        handle.write("\n")

    print(f"Analysis complete. Outputs written to: {output_dir}")
    print(f"Samples included: {len(samples)}")
    print("Network group sizes:")
    print(group_sizes.to_string())
    print("\nPCA first four PCs:")
    print(format_table(pca_variance.head(4)))
    print("\nPERMANOVA global:")
    print(format_table(global_permanova))
    print("\nPERMDISP global:")
    print(format_table(permdisp_global))
    print("\nKruskal-Wallis:")
    print(format_table(kruskal))


def main() -> None:
    script_path = Path(__file__).resolve()
    repo_root = script_path.parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=repo_root / "data" / "analysis.csv")
    parser.add_argument("--output", type=Path, default=repo_root / "results" / RUN_LABEL)
    parser.add_argument("--permutations", type=int, default=DEFAULT_PERMUTATIONS)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args()
    run(args.input.resolve(), args.output.resolve(), args.permutations, args.seed)


if __name__ == "__main__":
    main()
