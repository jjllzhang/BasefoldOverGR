#!/usr/bin/env python3
"""Plot benchmark metrics vs number of constraints from one or more benchmark CSV files."""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
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
from matplotlib.ticker import FuncFormatter, LogLocator


@dataclass(frozen=True)
class MetricSpec:
    key: str
    y_label: str
    output_tag: str


@dataclass
class SeriesData:
    label: str
    points_by_metric: Dict[str, List[Tuple[float, float]]]


@dataclass(frozen=True)
class CsvInputSpec:
    csv_path: Path
    include_families: Tuple[str, ...] = ()
    exclude_families: Tuple[str, ...] = ()
    include_context_ids: Tuple[str, ...] = ()
    exclude_context_ids: Tuple[str, ...] = ()


@dataclass(frozen=True)
class LegendEntry:
    label: str
    color: str
    marker: str


@dataclass(frozen=True)
class PlotResult:
    plot_path: Path
    legend_entries: Tuple[LegendEntry, ...]


@dataclass(frozen=True)
class FigureJob:
    output_prefix: str
    input_specs: Tuple[CsvInputSpec, ...]
    legend_mode: str = "auto"


METRIC_CATALOG = {
    "commit": MetricSpec(
        key="commit",
        y_label="Commit time (ms)",
        output_tag="commit_time_vs_number_of_constraints",
    ),
    "open": MetricSpec(
        key="open",
        y_label="Open time (ms)",
        output_tag="open_time_vs_number_of_constraints",
    ),
    "prover": MetricSpec(
        key="prover",
        y_label="Prover time (ms)",
        output_tag="prover_time_vs_number_of_constraints",
    ),
    "verifier": MetricSpec(
        key="verifier",
        y_label="Verifier time (ms)",
        output_tag="verifier_time_vs_number_of_constraints",
    ),
}

PROOF_SIZE_METRIC = MetricSpec(
    key="proof_size",
    y_label="Proof size (KB)",
    output_tag="proof_size_vs_number_of_constraints",
)

BACKEND_METRIC_COLUMNS = {
    "commit": "commit_mean_ms",
    "open": "open_mean_ms",
    "prover": "prove_mean_ms",
    "verifier": "verifier_mean_ms",
    "proof_size": "proof_size_kb",
}

COMPILER_METRIC_COLUMNS = {
    "commit": "commit_total_mean_ms",
    "open": "open_total_mean_ms",
    "prover": "prove_total_mean_ms",
    "verifier": "verify_total_mean_ms",
    "proof_size": "total_proof_size_kb",
}

DEFAULT_METRIC_NAMES = ["commit", "prover", "verifier", "proof_size"]


def _detect_metric_columns(fieldnames: Iterable[str]) -> Dict[str, str]:
    field_set = set(fieldnames)
    if set(BACKEND_METRIC_COLUMNS.values()).issubset(field_set):
        return BACKEND_METRIC_COLUMNS
    if set(COMPILER_METRIC_COLUMNS.values()).issubset(field_set):
        return COMPILER_METRIC_COLUMNS
    raise ValueError(
        "Unsupported CSV schema. Expected latest backend_eval_results.csv or "
        "compiler_eval_results.csv columns."
    )


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


def _parse_positive_int(raw: str | None) -> int | None:
    if raw is None:
        return None
    raw = raw.strip()
    if not raw or raw == "-":
        return None
    if not re.fullmatch(r"[1-9][0-9]*", raw):
        return None
    return int(raw)


def _parse_power_of_two(raw: str | None) -> int | None:
    value = _parse_positive_int(raw)
    if value is None:
        return None
    if value & (value - 1):
        raise ValueError(f"poly_dim must be a power of two, got {value}")
    return value


def _row_matches_spec(row: dict[str, str], input_spec: CsvInputSpec) -> bool:
    family = (row.get("family") or "").strip()
    context_id = (row.get("context_id") or "").strip()

    if input_spec.include_families and family not in input_spec.include_families:
        return False
    if input_spec.exclude_families and family in input_spec.exclude_families:
        return False
    if input_spec.include_context_ids and context_id not in input_spec.include_context_ids:
        return False
    if input_spec.exclude_context_ids and context_id in input_spec.exclude_context_ids:
        return False
    return True


_EXT_CHALLENGE_MARKER_RE = re.compile(
    r"\s*\(\s*ext(?:ension)?[- ]challenge\s*\)\s*",
    flags=re.IGNORECASE,
)
_GR_NOTATION_RE = re.compile(r"GR\(([^()]+)\)")
_GR_PAIR_RE = re.compile(r"GR\(\s*([^,;()]+)\s*[,;]\s*([^,;()]+)\s*\)")

_SPECIAL_FAMILY_LABELS = {
    "fri_based": "FRI-PCS",
    "ligero_based": "Ligero-PCS",
}

_FAMILY_DISPLAY_NAMES = {
    "frobenius": "Frobenius map",
    "ring_switch": "Ring-switching",
}

_BACKEND_CONTEXT_DOMAIN_LABELS = {
    "field-255": r"$\mathbb{F}_{2^{255}-19}$",
    "field-f2p256": r"$\mathbb{F}_{2^{256}}$",
    "field-prime64-ext": r"$\mathbb{F}_{2^{64}-59}$",
    "field-f2p64-ext": r"$\mathbb{F}_{2^{64}}$",
    "field-prime128-ext": r"$\mathbb{F}_{2^{128}-159}$",
    "field-f2p128-ext": r"$\mathbb{F}_{2^{128}}$",
    "field-f3p40-ext": r"$\mathbb{F}_{3^{40}}$",
    "field-f3p81-ext": r"$\mathbb{F}_{3^{81}}$",
    "ring-gr-2p16-162": r"$\mathrm{GR}(2^{16}, 162)$",
    "ring-gr-2p2-162": r"$\mathrm{GR}(2^{2}, 162)$",
    "ring-gr-2p16-64-ext": r"$\mathrm{GR}(2^{16}, 64)$",
    "ring-gr-2p16-128-ext": r"$\mathrm{GR}(2^{16}, 128)$",
    "ring-gr-2p2-64-ext": r"$\mathrm{GR}(2^{2}, 64)$",
    "ring-gr-2p2-128-ext": r"$\mathrm{GR}(2^{2}, 128)$",
}


def _strip_ext_challenge_marker(label: str) -> str:
    return _EXT_CHALLENGE_MARKER_RE.sub("", label).strip()


def _normalize_gr_notation(label: str) -> str:
    def _replace(match: re.Match[str]) -> str:
        inner = match.group(1)
        if ";" not in inner:
            return match.group(0)
        parts = [part.strip() for part in inner.split(";")]
        return f"GR({','.join(parts)})"

    return _GR_NOTATION_RE.sub(_replace, label)


def _context_series_label(csv_path: Path, labels: set[str]) -> str:
    if len(labels) == 1:
        raw_label = next(iter(labels))
    else:
        raw_label = csv_path.stem

    raw_label = _strip_ext_challenge_marker(raw_label)
    return _normalize_gr_notation(raw_label)


def _display_family_name(raw_family: str | None) -> str | None:
    if not raw_family:
        return None
    mapped_name = _FAMILY_DISPLAY_NAMES.get(raw_family)
    if mapped_name is not None:
        return mapped_name
    words = raw_family.replace("_", " ").split()
    if not words:
        return None
    return " ".join(word.capitalize() for word in words)


def _math_gr_label(raw_label: str) -> str | None:
    normalized = _normalize_gr_notation(_strip_ext_challenge_marker(raw_label))
    match = _GR_PAIR_RE.search(normalized)
    if match is None:
        return None

    left = re.sub(r"\^(\d+)", r"^{\1}", match.group(1).strip())
    right = match.group(2).strip()
    return rf"$\mathrm{{GR}}({left}, {right})$"


def _math_z2k_label(raw_label: str) -> str | None:
    normalized = _normalize_gr_notation(_strip_ext_challenge_marker(raw_label))
    match = _GR_PAIR_RE.search(normalized)
    if match is None:
        return None

    left = re.sub(r"\^(\d+)", r"^{\1}", match.group(1).strip())
    return rf"$\mathbb{{Z}}_{{{left}}}$"


def _stable_values(rows: Sequence[dict[str, str]], field: str) -> set[str]:
    return {
        (row.get(field) or "").strip()
        for row in rows
        if (row.get(field) or "").strip()
    }


def _format_split_field_value(csv_path: Path, field: str, value: str) -> str:
    if field == "display_name":
        return value
    if field == "family":
        return _display_family_name(value) or value
    if field in ("context_label", "context_id"):
        return _context_series_label(csv_path, {value})
    if field == "mode":
        return value
    return value


def _backend_series_label(csv_path: Path, rows: Sequence[dict[str, str]]) -> str | None:
    if _stable_values(rows, "family") != {"basefold"}:
        return None

    context_ids = _stable_values(rows, "context_id")
    if len(context_ids) != 1:
        return None
    context_id = next(iter(context_ids))

    domain_label = _BACKEND_CONTEXT_DOMAIN_LABELS.get(context_id)
    if domain_label is not None:
        return f"Basefold over {domain_label}"

    context_labels = _stable_values(rows, "context_label")
    if len(context_labels) == 1:
        return f"Basefold over {_context_series_label(csv_path, context_labels)}"
    return None


def _compiler_series_label(csv_path: Path, rows: Sequence[dict[str, str]]) -> str | None:
    families = _stable_values(rows, "family")
    if len(families) != 1:
        return None

    family = next(iter(families))
    if family not in {"ring_switch", "frobenius"}:
        return None

    family_label = _display_family_name(family)
    if not family_label:
        return None

    kappas = _stable_values(rows, "compiler_kappa")
    kappa_suffix = ""
    if len(kappas) == 1:
        kappa_suffix = rf" ($\kappa$ = {next(iter(kappas))})"

    context_labels = _stable_values(rows, "context_label")
    if len(context_labels) == 1:
        context_label = next(iter(context_labels))
        z2k_label = _math_z2k_label(context_label)
        if z2k_label is not None:
            return f"{family_label} over {z2k_label}{kappa_suffix}"
        return f"{family_label} over {_context_series_label(csv_path, context_labels)}{kappa_suffix}"

    context_ids = _stable_values(rows, "context_id")
    if len(context_ids) == 1:
        return f"{family_label} over {_context_series_label(csv_path, context_ids)}{kappa_suffix}"
    return f"{family_label}{kappa_suffix}"


def _special_series_label(rows: Sequence[dict[str, str]]) -> str | None:
    families = _stable_values(rows, "family")
    if len(families) != 1:
        return None

    family = next(iter(families))
    family_prefix = _SPECIAL_FAMILY_LABELS.get(family)
    if family_prefix is None:
        return None

    context_labels = _stable_values(rows, "context_label")
    if len(context_labels) != 1:
        return family_prefix

    gr_label = _math_gr_label(next(iter(context_labels)))
    if gr_label is None:
        return family_prefix
    return f"{family_prefix} over {gr_label}"


def _series_label(
    csv_path: Path,
    rows: Sequence[dict[str, str]],
    split_fields: Sequence[str],
    include_source_label: bool = False,
) -> str:
    special_label = _special_series_label(rows)
    if special_label is not None:
        label = special_label
    else:
        backend_label = _backend_series_label(csv_path, rows)
        if backend_label is not None:
            label = backend_label
        else:
            compiler_label = _compiler_series_label(csv_path, rows)
            if compiler_label is not None:
                label = compiler_label
            else:
                parts: List[str] = []
                for field in split_fields:
                    values = _stable_values(rows, field)
                    if len(values) != 1:
                        continue
                    rendered = _format_split_field_value(csv_path, field, next(iter(values)))
                    if rendered:
                        parts.append(rendered)

                if not parts:
                    display_names = _stable_values(rows, "display_name")
                    if len(display_names) == 1:
                        parts.append(next(iter(display_names)))
                    else:
                        families = _stable_values(rows, "family")
                        if len(families) == 1:
                            family_label = _display_family_name(next(iter(families)))
                            if family_label:
                                parts.append(family_label)

                label = " | ".join(parts) if parts else csv_path.stem

    if include_source_label:
        return f"{csv_path.stem} | {label}"
    return label


def _select_split_fields(rows: Sequence[dict[str, str]]) -> tuple[str, ...]:
    split_fields: List[str] = []
    for field in ("display_name", "family", "context_id", "mode"):
        values = _stable_values(rows, field)
        if len(values) > 1:
            split_fields.append(field)
    return tuple(split_fields)


def _format_power_of_two(value: float, _pos: int) -> str:
    if value <= 0:
        return ""
    exponent = math.log2(value)
    rounded = round(exponent)
    if abs(exponent - rounded) <= 1e-9:
        return rf"$2^{{{int(rounded)}}}$"
    return ""


def _format_power_of_ten(value: float, _pos: int) -> str:
    if value <= 0:
        return ""
    exponent = math.log10(value)
    rounded = round(exponent)
    if abs(exponent - rounded) <= 1e-9:
        return rf"$10^{{{int(rounded)}}}$"
    return ""


def load_series_from_csv(
    input_spec: CsvInputSpec,
    metric_keys: Iterable[str],
    proof_size_column: str | None = None,
    include_source_label: bool = False,
) -> List[SeriesData]:
    csv_path = input_spec.csv_path
    metric_keys = list(metric_keys)
    valid_rows: List[dict[str, str]] = []

    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"CSV has no header: {csv_path}")

        metric_columns = _detect_metric_columns(reader.fieldnames)
        resolved_columns = {
            metric_key: (
                proof_size_column
                if metric_key == "proof_size" and proof_size_column is not None
                else metric_columns[metric_key]
            )
            for metric_key in metric_keys
        }
        missing_columns = [
            column_name
            for column_name in resolved_columns.values()
            if column_name not in reader.fieldnames
        ]
        if missing_columns:
            raise ValueError(
                f"{csv_path}: missing expected columns: {', '.join(sorted(set(missing_columns)))}"
            )

        for row_index, row in enumerate(reader, start=2):
            status = (row.get("status") or "ok").strip().lower()
            if status not in ("ok", "success", ""):
                continue
            if not _row_matches_spec(row, input_spec):
                continue

            try:
                x_value = _parse_power_of_two(row.get("poly_dim"))
            except ValueError as exc:
                raise ValueError(
                    f"{csv_path}:{row_index}: invalid poly_dim: {exc}"
                ) from exc
            if x_value is None:
                continue

            normalized_row = dict(row)
            normalized_row["__x_value"] = str(x_value)
            valid_rows.append(normalized_row)

    split_fields = _select_split_fields(valid_rows)
    grouped_rows: Dict[tuple[str, ...], List[dict[str, str]]] = defaultdict(list)
    for row in valid_rows:
        group_key = tuple((row.get(field) or "").strip() for field in split_fields)
        grouped_rows[group_key].append(row)

    series_list: List[SeriesData] = []
    for group_key in sorted(grouped_rows):
        rows = grouped_rows[group_key]
        metric_values: Dict[str, Dict[float, List[float]]] = {
            key: defaultdict(list) for key in metric_keys
        }
        for row in rows:
            x_value = float(row["__x_value"])
            for metric_key in metric_keys:
                raw_value = row.get(resolved_columns[metric_key])
                metric_value = _parse_float(raw_value)
                if metric_value is not None:
                    metric_values[metric_key][x_value].append(metric_value)

        points_by_metric: Dict[str, List[Tuple[float, float]]] = {}
        for metric_key, by_x_value in metric_values.items():
            points = [
                (x_value, mean(values)) for x_value, values in by_x_value.items()
            ]
            points.sort(key=lambda pair: pair[0])
            points_by_metric[metric_key] = points

        series_list.append(
            SeriesData(
                label=_series_label(
                    csv_path,
                    rows,
                    split_fields,
                    include_source_label=include_source_label,
                ),
                points_by_metric=points_by_metric,
            )
        )

    return series_list


def plot_metric_vs_constraints(
    series_list: Sequence[SeriesData],
    metric: MetricSpec,
    output_path: Path,
    split_legend: bool = True,
    inline_legend_ncols: int | None = None,
    compiler_overlap_style: bool = False,
) -> PlotResult | None:
    fig, ax = plt.subplots(figsize=(10, 6))
    all_x_values: List[float] = []
    all_y_values: List[float] = []
    legend_entries: List[LegendEntry] = []
    visible_series = [
        (series, series.points_by_metric.get(metric.key, []))
        for series in series_list
        if series.points_by_metric.get(metric.key, [])
    ]
    plotted_count = len(visible_series)

    marker_cycle = cycle(["s", "o", "D", "v", "^", "P", "X", "<", ">"])

    for idx, (series, points) in enumerate(visible_series):
        x_values = [x for x, _ in points]
        y_values = [y for _, y in points]
        all_x_values.extend(x_values)
        all_y_values.extend(y for y in y_values if y > 0)

        marker = next(marker_cycle)
        linestyle = "-"
        linewidth = 2
        markersize = 6
        markerfacecolor = "white"
        markeredgewidth = 1.4
        alpha = 0.95
        plot_x_values = x_values
        zorder = 3 + idx
        if compiler_overlap_style:
            # Keep coordinates exact while separating nearly identical
            # compiler curves via line style and marker treatment.
            lower_label = series.label.lower()
            if "ring switch" in lower_label or "ring-switch" in lower_label:
                linestyle = (0, (4.5, 2.2))
                linewidth = 2.2
                markersize = 6.8
                markerfacecolor = "none"
                markeredgewidth = 1.8
                alpha = 0.98
                zorder = 10 + idx
            else:
                linestyle = "-"
                linewidth = 2.4
                markersize = 7.2
                markerfacecolor = "white"
                markeredgewidth = 1.6
                alpha = 0.9
        (line,) = ax.plot(
            plot_x_values,
            y_values,
            marker=marker,
            linewidth=linewidth,
            markersize=markersize,
            markerfacecolor=markerfacecolor,
            markeredgewidth=markeredgewidth,
            linestyle=linestyle,
            alpha=alpha,
            label=series.label,
            zorder=zorder,
        )
        line.set_path_effects(
            [pe.Stroke(linewidth=3.2, foreground="white", alpha=0.85), pe.Normal()]
        )
        legend_entries.append(
            LegendEntry(
                label=series.label,
                color=str(line.get_color()),
                marker=marker,
            )
        )
        plotted_count += 1

    if plotted_count == 0:
        plt.close(fig)
        return None

    ax.set_xlabel("Number of constraints")
    ax.set_ylabel(metric.y_label)
    ax.grid(True, linestyle="--", alpha=0.4)

    if all_x_values:
        x_min = min(all_x_values)
        x_max = max(all_x_values)
        x_min_pow = 2 ** math.floor(math.log2(x_min))
        x_max_pow = 2 ** math.ceil(math.log2(x_max))
        if compiler_overlap_style:
            x_min_pow = max(x_min_pow, 2 ** 7)
        if x_max_pow <= x_min_pow:
            x_max_pow = x_min_pow * 2
        ax.set_xlim(x_min_pow, x_max_pow)

        power_ticks = sorted(
            {
                int(round(x))
                for x in all_x_values
                if x > 0 and abs(math.log2(x) - round(math.log2(x))) <= 1e-8
            }
        )
        if power_ticks:
            power_ticks = [tick for tick in power_ticks if tick >= x_min_pow]
        if power_ticks:
            ax.set_xticks(power_ticks)
    ax.set_xscale("log", base=2)
    ax.xaxis.set_major_locator(LogLocator(base=2))
    ax.xaxis.set_major_formatter(FuncFormatter(_format_power_of_two))
    if all_y_values:
        y_min = min(all_y_values)
        y_max = max(all_y_values)
        y_min_pow = 10 ** math.floor(math.log10(y_min))
        y_max_pow = 10 ** math.ceil(math.log10(y_max))
        if y_max_pow <= y_min_pow:
            y_max_pow = y_min_pow * 10
        ax.set_ylim(y_min_pow, y_max_pow)
    ax.set_yscale("log", base=10)
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_major_formatter(FuncFormatter(_format_power_of_ten))

    if not split_legend and legend_entries:
        legend_cols = (
            inline_legend_ncols
            if inline_legend_ncols is not None
            else (1 if len(legend_entries) == 1 else 2 if len(legend_entries) <= 4 else 3)
        )
        ax.legend(
            handler_map={Line2D: HandlerLine2D(numpoints=1)},
            loc="upper left",
            ncol=legend_cols,
            handlelength=2.2,
            handletextpad=0.6,
            columnspacing=1.2,
            borderaxespad=0.3,
            labelspacing=0.4,
            fontsize=9,
            frameon=False,
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    if split_legend:
        return PlotResult(plot_path=output_path, legend_entries=tuple(legend_entries))
    return PlotResult(plot_path=output_path, legend_entries=tuple())


def render_legend_image(
    legend_entries: Sequence[LegendEntry],
    output_path: Path,
) -> None:
    legend_cols = 1 if len(legend_entries) == 1 else 2 if len(legend_entries) <= 4 else 3
    legend_rows = math.ceil(len(legend_entries) / legend_cols)
    legend_height = 0.9 + 0.38 * max(0, legend_rows - 1)

    handles = [
        Line2D(
            [],
            [],
            color=entry.color,
            marker=entry.marker,
            linewidth=2,
            markersize=6,
            markerfacecolor="white",
            markeredgewidth=1.4,
            linestyle="-",
            alpha=0.95,
        )
        for entry in legend_entries
    ]
    labels = [entry.label for entry in legend_entries]

    legend_fig = plt.figure(figsize=(10, legend_height))
    legend_fig.legend(
        handles,
        labels,
        handler_map={Line2D: HandlerLine2D(numpoints=1)},
        loc="center",
        ncol=legend_cols,
        handlelength=2.2,
        handletextpad=0.6,
        columnspacing=1.2,
        borderaxespad=0.0,
        labelspacing=0.4,
        fontsize=9,
        frameon=False,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    legend_fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(legend_fig)


def plot_benchmark_metrics(
    csv_paths: Sequence[Path] | None,
    output_dir: Path,
    output_prefix: str = "benchmark",
    proof_size_column: str | None = None,
    metrics_to_plot: Sequence[str] | None = None,
    legend_mode: str = "auto",
    input_specs: Sequence[CsvInputSpec] | None = None,
    include_source_label: bool | None = None,
) -> Tuple[List[Path], List[str]]:
    metric_catalog = {**METRIC_CATALOG, "proof_size": PROOF_SIZE_METRIC}
    selected_names = list(metrics_to_plot) if metrics_to_plot else DEFAULT_METRIC_NAMES
    metrics = [metric_catalog[name] for name in selected_names]
    metric_keys = [metric.key for metric in metrics]

    resolved_csv_paths = list(csv_paths or [])
    resolved_input_specs: Tuple[CsvInputSpec, ...]
    if input_specs is None:
        resolved_input_specs = tuple(CsvInputSpec(csv_path=path) for path in resolved_csv_paths)
    else:
        resolved_input_specs = tuple(input_specs)
        if not resolved_csv_paths:
            resolved_csv_paths = [spec.csv_path for spec in resolved_input_specs]

    if include_source_label is None:
        include_source_label = len(resolved_csv_paths) > 1

    split_legend = _should_split_legend(resolved_csv_paths, legend_mode)
    inline_legend_ncols = _inline_legend_ncols(resolved_csv_paths)
    compiler_overlap_style = _compiler_overlap_style(resolved_csv_paths)
    series_list: List[SeriesData] = []
    for input_spec in resolved_input_specs:
        series_list.extend(
            load_series_from_csv(
                input_spec,
                metric_keys,
                proof_size_column=proof_size_column,
                include_source_label=include_source_label,
            )
        )

    output_paths: List[Path] = []
    skipped_metrics: List[str] = []
    legend_cache: Dict[Tuple[LegendEntry, ...], Path] = {}
    stale_legend_candidates: List[Path] = []
    for metric in metrics:
        output_path = output_dir / f"{output_prefix}_{metric.output_tag}.png"
        legend_path = output_dir / f"{output_prefix}_{metric.output_tag}_legend.png"
        stale_legend_candidates.append(legend_path)

        plot_result = plot_metric_vs_constraints(
            series_list,
            metric,
            output_path,
            split_legend=split_legend,
            inline_legend_ncols=inline_legend_ncols,
            compiler_overlap_style=compiler_overlap_style,
        )
        if plot_result is not None:
            output_paths.append(plot_result.plot_path)

            if plot_result.legend_entries and plot_result.legend_entries not in legend_cache:
                render_legend_image(plot_result.legend_entries, legend_path)
                legend_cache[plot_result.legend_entries] = legend_path
                output_paths.append(legend_path)
        else:
            skipped_metrics.append(metric.key)

    keep_paths = set(output_paths)
    for stale_path in stale_legend_candidates:
        if stale_path not in keep_paths and stale_path.exists():
            stale_path.unlink()

    return output_paths, skipped_metrics


def _repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def _build_split_pcs_over_galois_ring_jobs() -> Tuple[FigureJob, ...]:
    results_new_dir = (_repo_root() / "results-new").resolve()
    fri_ligero_csv = (results_new_dir / "fri_ligero_based_eval_results.csv").resolve()
    jobs: List[FigureJob] = []
    for thread_dir, prefix_suffix in (("1-thread", "1_thread"), ("8-thread", "8_threads")):
        backend_csv = (results_new_dir / thread_dir / "backend_eval_results.csv").resolve()
        jobs.append(
            FigureJob(
                output_prefix=f"comparison_with_Basefold_{prefix_suffix}",
                input_specs=(
                    CsvInputSpec(
                        csv_path=backend_csv,
                        include_families=("basefold",),
                        exclude_context_ids=("ring-gr-2p16-162",),
                    ),
                ),
            )
        )
        jobs.append(
            FigureJob(
                output_prefix=f"comparison_with_JLXYY25_{prefix_suffix}",
                input_specs=(
                    CsvInputSpec(
                        csv_path=backend_csv,
                        include_families=("basefold",),
                        include_context_ids=("ring-gr-2p16-162",),
                    ),
                    CsvInputSpec(
                        csv_path=fri_ligero_csv,
                        include_families=("ligero_based",),
                    ),
                    CsvInputSpec(
                        csv_path=fri_ligero_csv,
                        include_families=("fri_based",),
                    ),
                ),
                legend_mode="inline",
            )
        )
    return tuple(jobs)


def _delete_matching_plots(output_dir: Path, patterns: Sequence[str]) -> List[Path]:
    removed_paths: List[Path] = []
    seen_paths: set[Path] = set()
    for pattern in patterns:
        for path in sorted(output_dir.glob(pattern)):
            if path in seen_paths:
                continue
            path.unlink()
            removed_paths.append(path)
            seen_paths.add(path)
    return removed_paths


def generate_split_pcs_over_galois_ring_figures(
    output_dir: Path,
    metrics_to_plot: Sequence[str] | None = None,
) -> Tuple[List[Path], List[str], List[Path]]:
    output_dir = output_dir.resolve()
    removed_paths = _delete_matching_plots(
        output_dir,
        (
            "PCSoverGaloisRing_*.png",
            "comparison_with_Basefold_*.png",
            "comparison_with_JLXYY25_*.png",
        ),
    )

    output_paths: List[Path] = []
    skipped_metrics: List[str] = []
    for job in _build_split_pcs_over_galois_ring_jobs():
        job_paths, job_skipped = plot_benchmark_metrics(
            csv_paths=[spec.csv_path for spec in job.input_specs],
            input_specs=job.input_specs,
            output_dir=output_dir,
            output_prefix=job.output_prefix,
            metrics_to_plot=metrics_to_plot,
            legend_mode=job.legend_mode,
            include_source_label=False,
        )
        output_paths.extend(job_paths)
        skipped_metrics.extend(
            f"{job.output_prefix}:{metric_name}" for metric_name in job_skipped
        )

    return output_paths, skipped_metrics, removed_paths


def _default_output_prefix(csv_paths: Sequence[Path], requested_prefix: str) -> str:
    if requested_prefix != "benchmark":
        return requested_prefix
    if len(csv_paths) != 1:
        return requested_prefix

    filename = csv_paths[0].name
    if filename == "backend_eval_results.csv":
        return "PCSoverGaloisRing"
    if filename == "compiler_eval_results.csv":
        return "PCSoverZ2K"
    return requested_prefix


def _should_split_legend(csv_paths: Sequence[Path], legend_mode: str) -> bool:
    if legend_mode == "split":
        return True
    if legend_mode == "inline":
        return False

    return not (
        len(csv_paths) == 1
        and csv_paths[0].name in {"backend_eval_results.csv", "compiler_eval_results.csv"}
    )


def _inline_legend_ncols(csv_paths: Sequence[Path]) -> int | None:
    if len(csv_paths) != 1:
        return None
    if csv_paths[0].name == "backend_eval_results.csv":
        return 3
    return None


def _compiler_overlap_style(csv_paths: Sequence[Path]) -> bool:
    return len(csv_paths) == 1 and csv_paths[0].name == "compiler_eval_results.csv"


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot commit/open/prover/verifier/proof-size vs number of constraints from one or more benchmark CSV files. Each line in a figure is one experiment, and a single CSV is automatically split when display_name, family, context_id, or mode varies."
        )
    )
    parser.add_argument("inputs", nargs="*", help="Input CSV file paths.")
    parser.add_argument(
        "-o",
        "--output-dir",
        default="result/plots",
        help="Output directory for PNG files (default: result/plots).",
    )
    parser.add_argument(
        "--prefix",
        default="benchmark",
        help=(
            "Prefix for output file names. If omitted, single-file backend and "
            "compiler runs default to PCSoverGaloisRing and PCSoverZ2K "
            "respectively; otherwise the default is benchmark."
        ),
    )
    parser.add_argument(
        "--proof-size-column",
        "--communication-column",
        dest="proof_size_column",
        default=None,
        help=(
            "Column to use for proof size "
            "(default: auto-select proof_size_kb for backend_eval_results.csv and total_proof_size_kb for compiler_eval_results.csv)."
        ),
    )
    parser.add_argument(
        "--metrics",
        nargs="+",
        help=(
            "Metrics to plot. Use names from: commit, open, prover, verifier, proof_size. "
            "Supports comma-separated tokens."
        ),
    )
    parser.add_argument(
        "--legend-mode",
        choices=("auto", "inline", "split"),
        default="auto",
        help=(
            "Legend placement mode. 'auto' keeps the previous behavior, "
            "'inline' always draws the legend inside the plot, and 'split' "
            "always writes a separate *_legend.png."
        ),
    )
    parser.add_argument(
        "--generate-split-pcs-over-galois-ring-figures",
        action="store_true",
        help=(
            "Regenerate the tracked results-new/figures comparison plots for "
            "PCS over Galois rings, split by thread count into Basefold and "
            "JLXYY25 comparison groups."
        ),
    )
    return parser.parse_args()


def _select_metrics(args: argparse.Namespace) -> List[str]:
    default_metrics = DEFAULT_METRIC_NAMES
    alias_map = {
        "commit": "commit",
        "commit_time": "commit",
        "open": "open",
        "open_time": "open",
        "prover": "prover",
        "prover_time": "prover",
        "verifier": "verifier",
        "verifier_time": "verifier",
        "proof_size": "proof_size",
        "proofsize": "proof_size",
        "proof": "proof_size",
        "communication": "proof_size",
    }

    if not args.metrics:
        return default_metrics

    tokens: List[str] = []
    for raw_item in args.metrics:
        tokens.extend(part.strip() for part in raw_item.split(",") if part.strip())

    selected: List[str] = []
    unknown: List[str] = []
    for token in tokens:
        normalized = token.lower().replace("-", "_")
        if normalized == "all":
            return default_metrics
        canonical = alias_map.get(normalized)
        if canonical is None:
            unknown.append(token)
            continue
        if canonical not in selected:
            selected.append(canonical)

    if unknown:
        valid = ", ".join(default_metrics + ["all"])
        raise ValueError(f"Unknown metric(s): {', '.join(unknown)}. Valid options: {valid}")
    if not selected:
        raise ValueError("No valid metrics selected.")
    return selected


def main() -> int:
    args = _parse_args()
    selected_metrics = _select_metrics(args)
    if args.generate_split_pcs_over_galois_ring_figures:
        if args.inputs:
            raise ValueError(
                "--generate-split-pcs-over-galois-ring-figures does not accept positional input CSVs."
            )

        default_output_dir = Path("result/plots")
        if Path(args.output_dir) == default_output_dir:
            output_dir = (_repo_root() / "results-new" / "figures").resolve()
        else:
            output_dir = Path(args.output_dir).expanduser().resolve()

        output_paths, skipped_metrics, removed_paths = generate_split_pcs_over_galois_ring_figures(
            output_dir=output_dir,
            metrics_to_plot=selected_metrics,
        )

        if removed_paths:
            print("Removed plots:")
            for removed_path in removed_paths:
                print(removed_path)

        if output_paths:
            print("Generated plots:")
            for output_path in output_paths:
                print(output_path)
        else:
            print("Generated plots: none")

        if skipped_metrics:
            print("Skipped metrics with no valid data:")
            for metric in skipped_metrics:
                print(metric)

        return 0

    if not args.inputs:
        raise ValueError("Provide at least one input CSV, or use --generate-split-pcs-over-galois-ring-figures.")

    csv_paths = [Path(path).expanduser().resolve() for path in args.inputs]
    missing_files = [str(path) for path in csv_paths if not path.exists()]
    if missing_files:
        raise FileNotFoundError(f"Input file(s) not found: {', '.join(missing_files)}")

    output_dir = Path(args.output_dir).expanduser().resolve()
    output_prefix = _default_output_prefix(csv_paths, args.prefix)
    output_paths, skipped_metrics = plot_benchmark_metrics(
        csv_paths=csv_paths,
        output_dir=output_dir,
        output_prefix=output_prefix,
        proof_size_column=args.proof_size_column,
        metrics_to_plot=selected_metrics,
        legend_mode=args.legend_mode,
    )

    if output_paths:
        print("Generated plots:")
        for output_path in output_paths:
            print(output_path)
    else:
        print("Generated plots: none")

    if skipped_metrics:
        print("Skipped metrics with no valid data:")
        for metric in skipped_metrics:
            print(metric)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
