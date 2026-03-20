from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


def test_reproduce_e001_script_generates_expected_outputs(tmp_path):
    repo_root = Path(__file__).resolve().parents[3]
    script_path = repo_root / "python" / "crossmapy" / "scripts" / "reproduce_e001.py"
    data_dir = repo_root / "Spatial convergent cross mapping" / "3561957"
    out_dir = tmp_path / "e001-results"

    result = subprocess.run(
        [
            sys.executable,
            str(script_path),
            "--data-dir",
            str(data_dir),
            "--iterations",
            "5",
            "--out-dir",
            str(out_dir),
            "--seed",
            "2718",
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr

    summary_csv = out_dir / "summary.csv"
    summary_json = out_dir / "summary.json"
    plot_low = out_dir / "rho_vs_L_lvl0.png"
    plot_high = out_dir / "rho_vs_L_lvl3.png"

    assert summary_csv.exists()
    assert summary_json.exists()
    assert plot_low.exists()
    assert plot_high.exists()

    payload = json.loads(summary_json.read_text())
    assert payload["iterations"] == 5
    assert "lvl0" in payload["scenarios"]
    assert "lvl3" in payload["scenarios"]
    assert "python" in payload["scenarios"]["lvl0"]
    assert "r" in payload["scenarios"]["lvl0"]
