#!/usr/bin/env python3
"""Plot benchmark metrics vs d from one or more benchmark CSV files."""

from __future__ import annotations

import argparse
import csv
import math
import os
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from itertools import cycle
from pathlib import Path
from statistics import mean
from typing import Dict, Iterable, List, Sequence, Tuple

# Avoid Matplotlib warnings when HOME config dir is not writable.
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.legend_handler import HandlerLine2D
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator


@dataclass(frozen=True)
class MetricSpec:
    key: str
    y_label: str
    output_tag: str


@dataclass
class SeriesData:
    label: str
    points_by_metric: Dict[str, List[Tuple[float, float]]]


def _parse_float(raw: str | None) -> float | None:
    if raw is None:
        return None
    raw = raw.strip()
    if not raw or raw == "-":
        return None
    try:
        value = float(raw)
    except ValueError:
        return None
    if math.isnan(value) or math.isinf(value):
        return None
    return value


def _series_label(csv_path: Path, labels: set[str]) -> str:
    if len(labels) == 1:
        return next(iter(labels))
    return csv_path.stem


def load_series_from_csv(csv_path: Path, metric_keys: Iterable[str]) -> SeriesData:
    metric_values: Dict[str, Dict[float, List[float]]] = {
        key: defaultdict(list) for key in metric_keys
    }
    labels: set[str] = set()

    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"CSV has no header: {csv_path}")

        for row in reader:
            status = (row.get("status") or "ok").strip().lower()
            if status not in ("ok", "success", ""):
                continue

            label = (row.get("context_label") or row.get("context_id") or "").strip()
            if label:
                labels.add(label)

            d = _parse_float(row.get("d"))
            if d is None:
                poly_dim = _parse_float(row.get("poly_dim"))
                if poly_dim is not None and poly_dim > 0:
                    d = math.log2(poly_dim)
            if d is None:
                continue

            for metric_key in metric_keys:
                metric_value = _parse_float(row.get(metric_key))
                if metric_value is not None:
                    metric_values[metric_key][d].append(metric_value)

    points_by_metric: Dict[str, List[Tuple[float, float]]] = {}
    for metric_key, by_d in metric_values.items():
        points = [(d, mean(values)) for d, values in by_d.items()]
        points.sort(key=lambda pair: pair[0])
        points_by_metric[metric_key] = points

    return SeriesData(label=_series_label(csv_path, labels), points_by_metric=points_by_metric)


def plot_metric_vs_d(
    series_list: Sequence[SeriesData], metric: MetricSpec, output_path: Path
) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    plotted_count = 0
    all_x_values: List[float] = []

    marker_cycle = cycle(["s", "o", "D", "v", "^", "P", "X", "<", ">"])

    for idx, series in enumerate(series_list):
        points = series.points_by_metric.get(metric.key, [])
        if not points:
            continue
        x_values = [x for x, _ in points]
        y_values = [y for _, y in points]
        all_x_values.extend(x_values)

        marker = next(marker_cycle)
        (line,) = ax.plot(
            x_values,
            y_values,
            marker=marker,
            linewidth=2,
            markersize=6,
            markerfacecolor="white",
            markeredgewidth=1.4,
            linestyle="-",
            alpha=0.95,
            markevery=1,
            label=series.label,
            zorder=3 + idx,
        )
        line.set_path_effects(
            [pe.Stroke(linewidth=3.2, foreground="white", alpha=0.85), pe.Normal()]
        )
        plotted_count += 1

    if plotted_count == 0:
        plt.close(fig)
        raise ValueError(f"No valid data for metric '{metric.key}'")

    ax.set_xlabel("log2(# of polynomial dimensions)")
    ax.set_ylabel(metric.y_label)
    ax.grid(True, linestyle="--", alpha=0.4)

    if all_x_values:
        integer_ticks = sorted(
            {
                int(round(x))
                for x in all_x_values
                if abs(x - round(x)) <= 1e-8
            }
        )
        if integer_ticks:
            ax.set_xticks(integer_ticks)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    if plotted_count > 1:
        ax.legend(
            handler_map={Line2D: HandlerLine2D(numpoints=1)},
            loc="upper left",
            handlelength=2.2,
            handletextpad=0.6,
            borderaxespad=0.3,
            labelspacing=0.4,
        )

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_benchmark_metrics(
    csv_paths: Sequence[Path],
    output_dir: Path,
    output_prefix: str = "benchmark",
    communication_column: str = "proof_size_kb",
) -> List[Path]:
    metrics = [
        MetricSpec(
            key="commit_mean_ms",
            y_label="Commit time (ms)",
            output_tag="commit_time_vs_d",
        ),
        MetricSpec(
            key="prover_mean_ms",
            y_label="Prover time (ms)",
            output_tag="prover_time_vs_d",
        ),
        MetricSpec(
            key="verifier_mean_ms",
            y_label="Verifier time (ms)",
            output_tag="verifier_time_vs_d",
        ),
        MetricSpec(
            key=communication_column,
            y_label="communication (KB)",
            output_tag="communication_vs_d",
        ),
    ]
    metric_keys = [metric.key for metric in metrics]
    series_list = [load_series_from_csv(csv_path, metric_keys) for csv_path in csv_paths]

    output_paths: List[Path] = []
    for metric in metrics:
        output_path = output_dir / f"{output_prefix}_{metric.output_tag}.png"
        plot_metric_vs_d(series_list, metric, output_path)
        output_paths.append(output_path)

    return output_paths


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot commit/prover/verifier/communication vs d from one or more benchmark CSV files."
        )
    )
    parser.add_argument("inputs", nargs="+", help="Input CSV file paths.")
    parser.add_argument(
        "-o",
        "--output-dir",
        default="result/plots",
        help="Output directory for PNG files (default: result/plots).",
    )
    parser.add_argument(
        "--prefix",
        default="benchmark",
        help="Prefix for output file names (default: benchmark).",
    )
    parser.add_argument(
        "--communication-column",
        "--proof-size-column",
        dest="communication_column",
        default="proof_size_kb",
        help="Column to use for communication (default: proof_size_kb).",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    csv_paths = [Path(path).expanduser().resolve() for path in args.inputs]

    missing_files = [str(path) for path in csv_paths if not path.exists()]
    if missing_files:
        raise FileNotFoundError(f"Input file(s) not found: {', '.join(missing_files)}")

    output_dir = Path(args.output_dir).expanduser().resolve()
    output_paths = plot_benchmark_metrics(
        csv_paths=csv_paths,
        output_dir=output_dir,
        output_prefix=args.prefix,
        communication_column=args.communication_column,
    )

    print("Generated plots:")
    for output_path in output_paths:
        print(output_path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
