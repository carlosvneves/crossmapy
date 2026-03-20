#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[3]
PACKAGE_SRC = SCRIPT_PATH.parents[1] / "src"
if str(PACKAGE_SRC) not in sys.path:
    sys.path.insert(0, str(PACKAGE_SRC))

from crossmapy import CCM_boot, SSR_pred_boot, ccmtest  # noqa: E402


ARTICLE_TARGETS = {
    "lvl0": {"A_to_B": 0.614, "B_to_A": 0.004},
    "lvl3": {"A_to_B": 0.221, "B_to_A": 0.643},
}


@dataclass
class ScenarioData:
    name: str
    A: np.ndarray
    B: np.ndarray


def _find_header_row_index(path: Path) -> int:
    with path.open("r", encoding="utf-8") as handle:
        for idx, line in enumerate(handle):
            line_norm = line.strip().strip('"').lower()
            if line_norm.startswith("index") and "agropyron" in line_norm:
                return idx
    raise ValueError(f"Could not find CSV header in {path}")


def _to_float_or_nan(value: str) -> float:
    value = (value or "").strip()
    if not value:
        return float("nan")
    if value.upper() == "NA":
        return float("nan")
    return float(value)


def load_e001_series(path: Path, scenario_name: str) -> ScenarioData:
    header_idx = _find_header_row_index(path)
    with path.open("r", encoding="utf-8") as handle:
        for _ in range(header_idx):
            next(handle)
        reader = csv.DictReader(handle)
        A_vals = []
        B_vals = []
        for row in reader:
            A_vals.append(_to_float_or_nan(row.get("Agropyron repens", "")))
            B_vals.append(_to_float_or_nan(row.get("Schizachyrium scoparium", "")))
    return ScenarioData(name=scenario_name, A=np.array(A_vals), B=np.array(B_vals))


def choose_embedding_dimension(series: np.ndarray, max_e: int = 16) -> int:
    best_e = None
    best_rho = -np.inf
    for e_val in range(2, max_e + 1):
        out = SSR_pred_boot(A=series, E=e_val, tau=1, predstep=1)
        rho = float(out["rho"])
        if np.isfinite(rho) and rho > best_rho:
            best_rho = rho
            best_e = e_val
    if best_e is None:
        raise ValueError("Could not find finite rho when selecting embedding dimension.")
    return best_e


def _serialize_series_to_csv(path: Path, A: np.ndarray, B: np.ndarray) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["A", "B"])
        for av, bv in zip(A, B):
            a_txt = "" if np.isnan(av) else f"{av:.17g}"
            b_txt = "" if np.isnan(bv) else f"{bv:.17g}"
            writer.writerow([a_txt, b_txt])


def _parse_r_list(text: str) -> list[float]:
    if not text:
        return []
    out = []
    for token in text.split(","):
        tok = token.strip()
        if not tok or tok.upper() == "NA":
            out.append(float("nan"))
        else:
            out.append(float(tok))
    return out


def run_r_baseline(
    A: np.ndarray, B: np.ndarray, e_a: int, e_b: int, iterations: int
) -> dict[str, object]:
    with tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False, encoding="utf-8") as tmp:
        tmp_path = Path(tmp.name)
    try:
        _serialize_series_to_csv(tmp_path, A, B)
        r_code = r"""
args <- commandArgs(trailingOnly=TRUE)
csv_path <- args[[1]]
e_a <- as.integer(args[[2]])
e_b <- as.integer(args[[3]])
iterations <- as.integer(args[[4]])

suppressPackageStartupMessages(library(multispatialCCM))
dat <- read.csv(csv_path)
A <- dat$A
B <- dat$B

ab <- CCM_boot(A, B, E=e_a, tau=1, iterations=iterations)
ba <- CCM_boot(B, A, E=e_b, tau=1, iterations=iterations)
pv <- ccmtest(ab, ba)
pv_a <- if (is.list(pv)) pv$pval_a_cause_b else unname(pv["pval_a_cause_b"])
pv_b <- if (is.list(pv)) pv$pval_b_cause_a else unname(pv["pval_b_cause_a"])

cat("P_A_TO_B=", pv_a, "\n", sep="")
cat("P_B_TO_A=", pv_b, "\n", sep="")
cat("L_AB=", paste(ab$Lobs, collapse=","), "\n", sep="")
cat("RHO_AB=", paste(ab$rho, collapse=","), "\n", sep="")
cat("L_BA=", paste(ba$Lobs, collapse=","), "\n", sep="")
cat("RHO_BA=", paste(ba$rho, collapse=","), "\n", sep="")
"""
        proc = subprocess.run(
            ["Rscript", "-e", r_code, str(tmp_path), str(e_a), str(e_b), str(iterations)],
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(f"R baseline failed:\n{proc.stderr}")
        parsed: dict[str, str] = {}
        for line in proc.stdout.splitlines():
            if "=" in line:
                key, value = line.split("=", 1)
                parsed[key.strip()] = value.strip()
        return {
            "p_a_to_b": float(parsed["P_A_TO_B"]),
            "p_b_to_a": float(parsed["P_B_TO_A"]),
            "L_ab": np.array(_parse_r_list(parsed.get("L_AB", "")), dtype=float),
            "rho_ab": np.array(_parse_r_list(parsed.get("RHO_AB", "")), dtype=float),
            "L_ba": np.array(_parse_r_list(parsed.get("L_BA", "")), dtype=float),
            "rho_ba": np.array(_parse_r_list(parsed.get("RHO_BA", "")), dtype=float),
        }
    finally:
        tmp_path.unlink(missing_ok=True)


def run_python_ccm(
    A: np.ndarray, B: np.ndarray, e_a: int, e_b: int, iterations: int, seed: int
) -> dict[str, object]:
    np.random.seed(seed)
    ab = CCM_boot(A=A, B=B, E=e_a, tau=1, iterations=iterations)
    np.random.seed(seed)
    ba = CCM_boot(A=B, B=A, E=e_b, tau=1, iterations=iterations)
    pvals = ccmtest(ab, ba)
    return {
        "p_a_to_b": float(pvals["pval_a_cause_b"]),
        "p_b_to_a": float(pvals["pval_b_cause_a"]),
        "L_ab": np.asarray(ab["Lobs"], dtype=float),
        "rho_ab": np.asarray(ab["rho"], dtype=float),
        "L_ba": np.asarray(ba["Lobs"], dtype=float),
        "rho_ba": np.asarray(ba["rho"], dtype=float),
    }


def _rho_mae_on_common_l(
    L_python: np.ndarray, rho_python: np.ndarray, L_r: np.ndarray, rho_r: np.ndarray
) -> tuple[float, int]:
    py_map = {int(l): float(r) for l, r in zip(L_python, rho_python)}
    r_map = {int(l): float(r) for l, r in zip(L_r, rho_r)}
    common = sorted(set(py_map).intersection(r_map))
    if not common:
        return float("nan"), 0
    abs_errors = [abs(py_map[l] - r_map[l]) for l in common]
    return float(np.mean(abs_errors)), len(common)


def _json_default(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        v = float(value)
        if math.isnan(v):
            return None
        return v
    raise TypeError(f"Object not JSON serializable: {type(value)!r}")


def _plot_scenario(
    out_path: Path, scenario_key: str, py_out: dict[str, object], r_out: dict[str, object]
) -> None:
    plt.figure(figsize=(9, 5))
    plt.plot(py_out["L_ab"], py_out["rho_ab"], label="Python A->B", linewidth=2)
    plt.plot(py_out["L_ba"], py_out["rho_ba"], label="Python B->A", linewidth=2)
    plt.plot(r_out["L_ab"], r_out["rho_ab"], "--", label="R A->B", linewidth=2)
    plt.plot(r_out["L_ba"], r_out["rho_ba"], "--", label="R B->A", linewidth=2)
    plt.xlabel("Library Length (L)")
    plt.ylabel("Pearson Correlation (rho)")
    plt.title(f"E001 {scenario_key}: rho vs L")
    plt.legend(loc="best")
    plt.grid(alpha=0.2)
    text = (
        f"Article p(A->B)={ARTICLE_TARGETS[scenario_key]['A_to_B']:.3f}, "
        f"p(B->A)={ARTICLE_TARGETS[scenario_key]['B_to_A']:.3f}\n"
        f"Python p(A->B)={py_out['p_a_to_b']:.3f}, p(B->A)={py_out['p_b_to_a']:.3f}\n"
        f"R p(A->B)={r_out['p_a_to_b']:.3f}, p(B->A)={r_out['p_b_to_a']:.3f}"
    )
    plt.figtext(0.02, 0.01, text, ha="left", fontsize=9)
    plt.tight_layout(rect=(0, 0.08, 1, 1))
    plt.savefig(out_path, dpi=150)
    plt.close()


def _write_summary_csv(rows: list[dict[str, object]], out_path: Path) -> None:
    columns = [
        "scenario",
        "direction",
        "E_source",
        "E_response",
        "python_p_value",
        "r_p_value",
        "article_p_value",
        "abs_diff_python_r_p",
        "mae_rho_common_l",
        "n_common_l",
    ]
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reproduce E001 (nitrogen addition) and compare Python vs R baselines."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=REPO_ROOT / "Spatial convergent cross mapping" / "3561957",
        help="Directory containing E001 CSV files.",
    )
    parser.add_argument(
        "--iterations", type=int, default=1000, help="Bootstrap iterations for CCM."
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=REPO_ROOT / "python" / "crossmapy" / "results" / "e001",
        help="Output directory for summary files and plots.",
    )
    parser.add_argument("--seed", type=int, default=2718, help="Random seed.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    scenarios = {
        "lvl0": load_e001_series(args.data_dir / "e001_arssnlvl0.csv", "lvl0"),
        "lvl3": load_e001_series(args.data_dir / "e001_arssnlvl3.csv", "lvl3"),
    }

    summary_rows: list[dict[str, object]] = []
    json_payload: dict[str, object] = {
        "iterations": args.iterations,
        "seed": args.seed,
        "data_dir": str(args.data_dir),
        "scenarios": {},
    }

    for scenario_key, scenario in scenarios.items():
        e_a = choose_embedding_dimension(scenario.A)
        e_b = choose_embedding_dimension(scenario.B)

        py_out = run_python_ccm(
            A=scenario.A,
            B=scenario.B,
            e_a=e_a,
            e_b=e_b,
            iterations=args.iterations,
            seed=args.seed,
        )
        r_out = run_r_baseline(
            A=scenario.A, B=scenario.B, e_a=e_a, e_b=e_b, iterations=args.iterations
        )

        mae_ab, n_ab = _rho_mae_on_common_l(
            py_out["L_ab"], py_out["rho_ab"], r_out["L_ab"], r_out["rho_ab"]
        )
        mae_ba, n_ba = _rho_mae_on_common_l(
            py_out["L_ba"], py_out["rho_ba"], r_out["L_ba"], r_out["rho_ba"]
        )

        summary_rows.append(
            {
                "scenario": scenario_key,
                "direction": "A_to_B",
                "E_source": e_a,
                "E_response": e_b,
                "python_p_value": py_out["p_a_to_b"],
                "r_p_value": r_out["p_a_to_b"],
                "article_p_value": ARTICLE_TARGETS[scenario_key]["A_to_B"],
                "abs_diff_python_r_p": abs(py_out["p_a_to_b"] - r_out["p_a_to_b"]),
                "mae_rho_common_l": mae_ab,
                "n_common_l": n_ab,
            }
        )
        summary_rows.append(
            {
                "scenario": scenario_key,
                "direction": "B_to_A",
                "E_source": e_b,
                "E_response": e_a,
                "python_p_value": py_out["p_b_to_a"],
                "r_p_value": r_out["p_b_to_a"],
                "article_p_value": ARTICLE_TARGETS[scenario_key]["B_to_A"],
                "abs_diff_python_r_p": abs(py_out["p_b_to_a"] - r_out["p_b_to_a"]),
                "mae_rho_common_l": mae_ba,
                "n_common_l": n_ba,
            }
        )

        _plot_scenario(args.out_dir / f"rho_vs_L_{scenario_key}.png", scenario_key, py_out, r_out)

        json_payload["scenarios"][scenario_key] = {
            "E_A": e_a,
            "E_B": e_b,
            "article_targets": ARTICLE_TARGETS[scenario_key],
            "python": py_out,
            "r": r_out,
            "comparison": {
                "A_to_B": {"mae_rho_common_l": mae_ab, "n_common_l": n_ab},
                "B_to_A": {"mae_rho_common_l": mae_ba, "n_common_l": n_ba},
            },
        }

    _write_summary_csv(summary_rows, args.out_dir / "summary.csv")
    (args.out_dir / "summary.json").write_text(
        json.dumps(json_payload, indent=2, default=_json_default), encoding="utf-8"
    )

    print(f"Results written to: {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
