from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
import statistics
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
REFERENCE_ROOT = REPO_ROOT.parent / "reference_implementation"
REFERENCE_SCRIPTS = REPO_ROOT / "test_data" / "scripts"
DEFAULT_WORKER = REPO_ROOT / "build" / "Release" / "long_benchmark_worker.exe"
MODES = ("fast", "hybrid", "complete")
TEMPERATURES = ("cold", "warm")


@dataclass
class BenchmarkRow:
  sequence: int
  timestamp: str
  implementation: str
  temperature: str
  mode: str
  iteration: int
  construct_seconds: float
  solve_seconds: float
  total_seconds: float
  fixture: str
  note: str = ""


def finite_seconds(value: float) -> float:
  if not math.isfinite(value) or value < 0.0:
    raise RuntimeError(f"invalid timing value: {value}")
  return float(value)


def add_reference_imports() -> Any:
  sys.path.insert(0, str(REFERENCE_SCRIPTS))
  import generate_reference_fixtures as fixtures

  fixtures.add_reference_to_path()
  return fixtures


def load_reference_fixture(fixtures: Any, name: str) -> dict[str, Any]:
  for fixture in fixtures.fixture_definitions():
    if fixture["name"] == name:
      return fixture
  raise RuntimeError(f"unknown fixture: {name}")


def make_reference_solver(fixtures: Any, geodrape: Any, fixture: dict[str, Any], mode: str) -> Any:
  return geodrape.GeoDrapeSolver(
      fixtures.pyvista_from_vf(fixture["vertices"], fixture["faces"]),
      solver_options=fixture["solver_options"],
      default_fast=(mode == "fast"),
  )


def solve_reference_mode(solver: Any, fixture: dict[str, Any], mode: str) -> str:
  angle = fixture["angle_degrees"]
  seed_xy = fixture["seed_xy"]
  if mode == "fast":
    solver.solve(angle=angle, seed_pnt=seed_xy, fast=True)
    return ""
  if mode == "complete":
    solver.solve(angle=angle, seed_pnt=seed_xy, fast=False)
    return ""
  if mode == "hybrid":
    origin_bary, solver.dir_vecs = solver.place_origin(seed_xy, angle)
    source_curves = solver.trace_generators(origin_bary, solver.dir_vecs)
    fields = solver.solve_heat(source_curves)
    vectors = solver.solve_heat_fast(source_curves)
    solver.surf["dist_0"] = fields[0]
    solver.surf["dist_1"] = fields[1]
    solver.add_fields_fast(vectors)
    return "reference_hybrid_emulated_compute_distance_plus_compute_Y_vertex"
  raise RuntimeError(f"unknown reference mode: {mode}")


def run_reference_cold(fixtures: Any,
                       geodrape: Any,
                       fixture: dict[str, Any],
                       mode: str) -> tuple[float, float, float, str]:
  construct_start = time.perf_counter()
  solver = make_reference_solver(fixtures, geodrape, fixture, mode)
  construct_seconds = finite_seconds(time.perf_counter() - construct_start)

  solve_start = time.perf_counter()
  note = solve_reference_mode(solver, fixture, mode)
  solve_seconds = finite_seconds(time.perf_counter() - solve_start)
  return construct_seconds, solve_seconds, construct_seconds + solve_seconds, note


def prepare_reference_warm(fixtures: Any,
                           geodrape: Any,
                           fixture: dict[str, Any]) -> dict[str, Any]:
  solvers: dict[str, Any] = {}
  for mode in MODES:
    solver = make_reference_solver(fixtures, geodrape, fixture, mode)
    solve_reference_mode(solver, fixture, mode)
    solvers[mode] = solver
  return solvers


def run_reference_warm(solvers: dict[str, Any],
                       fixture: dict[str, Any],
                       mode: str) -> tuple[float, float, float, str]:
  solve_start = time.perf_counter()
  note = solve_reference_mode(solvers[mode], fixture, mode)
  solve_seconds = finite_seconds(time.perf_counter() - solve_start)
  return 0.0, solve_seconds, solve_seconds, note


class CppWorker:
  def __init__(self, executable: Path, fixture: str) -> None:
    if not executable.exists():
      raise RuntimeError(f"C++ worker executable not found: {executable}")
    self.process = subprocess.Popen(
        [str(executable), fixture],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    )
    ready = self._readline()
    if ready != "READY,worker":
      raise RuntimeError(f"unexpected C++ worker startup response: {ready}")

  def _readline(self) -> str:
    if self.process.stdout is None:
      raise RuntimeError("C++ worker stdout is closed")
    line = self.process.stdout.readline()
    if not line:
      stderr = self.process.stderr.read() if self.process.stderr else ""
      raise RuntimeError(f"C++ worker exited unexpectedly: {stderr}")
    return line.strip()

  def command(self, line: str) -> str:
    if self.process.stdin is None:
      raise RuntimeError("C++ worker stdin is closed")
    self.process.stdin.write(line + "\n")
    self.process.stdin.flush()
    return self._readline()

  def prepare_warm(self) -> None:
    for mode in MODES:
      response = self.command(f"PREPARE,warm,{mode}")
      if response != f"READY,{mode}":
        raise RuntimeError(f"unexpected C++ warm response: {response}")

  def run(self, temperature: str, mode: str, iteration: int) -> tuple[float, float, float, str]:
    response = self.command(f"RUN,{temperature},{mode},{iteration}")
    parts = response.split(",")
    if len(parts) != 8 or parts[0] != "RESULT":
      raise RuntimeError(f"unexpected C++ run response: {response}")
    return (
        finite_seconds(float(parts[5])),
        finite_seconds(float(parts[6])),
        finite_seconds(float(parts[7])),
        "",
    )

  def close(self) -> None:
    try:
      self.command("EXIT")
    finally:
      self.process.terminate()


def summarize(values: list[float]) -> dict[str, float]:
  ordered = sorted(values)
  return {
      "mean": statistics.fmean(ordered),
      "median": statistics.median(ordered),
      "min": ordered[0],
      "max": ordered[-1],
  }


def write_summary(rows: list[BenchmarkRow], path: Path) -> None:
  grouped: dict[tuple[str, str, str], list[float]] = {}
  for row in rows:
    grouped.setdefault((row.implementation, row.temperature, row.mode), []).append(row.total_seconds)

  summary = {
      "/".join(key): summarize(values)
      for key, values in sorted(grouped.items())
  }
  path.write_text(json.dumps(summary, indent=2), encoding="utf-8")


def write_plots(rows: list[BenchmarkRow], output_prefix: Path) -> None:
  try:
    import matplotlib.pyplot as plt
  except Exception as exc:
    print(f"plot skipped: matplotlib import failed: {exc}", file=sys.stderr)
    return

  for temperature in TEMPERATURES:
    plt.figure(figsize=(11, 7))
    for mode in MODES:
      for implementation in ("reference", "cpp"):
        series = [
            row.total_seconds
            for row in rows
            if row.temperature == temperature and
               row.mode == mode and
               row.implementation == implementation
        ]
        if not series:
          continue
        cumulative = []
        total = 0.0
        for value in series:
          total += value
          cumulative.append(total)
        plt.plot(
            range(1, len(cumulative) + 1),
            cumulative,
            label=f"{implementation} {mode}",
        )
    plt.xlabel("Measured run")
    plt.ylabel("Cumulative total seconds")
    plt.title(f"{temperature.capitalize()} solve cumulative time")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_prefix.with_name(output_prefix.name + f"_{temperature}.png"), dpi=160)
    plt.close()


def write_csv(rows: list[BenchmarkRow], path: Path) -> None:
  with path.open("w", newline="", encoding="utf-8") as file:
    writer = csv.DictWriter(file, fieldnames=list(BenchmarkRow.__dataclass_fields__.keys()))
    writer.writeheader()
    for row in rows:
      writer.writerow(row.__dict__)


def parse_args() -> argparse.Namespace:
  parser = argparse.ArgumentParser(description="Long alternating benchmark for Python reference and C++ intrinsic core.")
  parser.add_argument("--fixture", default="demo_part")
  parser.add_argument("--runs", type=int, default=12, help="Measured runs per implementation/temperature/mode.")
  parser.add_argument("--output-dir", type=Path, default=REPO_ROOT / "benchmark_runs")
  parser.add_argument("--worker", type=Path, default=DEFAULT_WORKER)
  parser.add_argument("--no-plot", action="store_true")
  return parser.parse_args()


def main() -> int:
  args = parse_args()
  if args.runs <= 0:
    raise RuntimeError("--runs must be positive")

  fixtures = add_reference_imports()
  import GeoDrape as geodrape

  fixture = load_reference_fixture(fixtures, args.fixture)
  args.output_dir.mkdir(parents=True, exist_ok=True)

  timestamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
  output_prefix = args.output_dir / f"long_benchmark_{args.fixture}_{timestamp}"
  csv_path = output_prefix.with_suffix(".csv")
  summary_path = output_prefix.with_suffix(".summary.json")

  worker = CppWorker(args.worker, args.fixture)
  rows: list[BenchmarkRow] = []
  sequence = 0

  try:
    reference_warm = prepare_reference_warm(fixtures, geodrape, fixture)
    worker.prepare_warm()

    for iteration in range(args.runs):
      for temperature_index, temperature in enumerate(TEMPERATURES):
        for mode_index, mode in enumerate(MODES):
          implementations = ["reference", "cpp"]
          if (iteration + temperature_index + mode_index) % 2:
            implementations.reverse()

          for implementation in implementations:
            if implementation == "reference":
              if temperature == "cold":
                construct, solve, total, note = run_reference_cold(fixtures, geodrape, fixture, mode)
              else:
                construct, solve, total, note = run_reference_warm(reference_warm, fixture, mode)
            else:
              construct, solve, total, note = worker.run(temperature, mode, iteration)

            rows.append(BenchmarkRow(
                sequence=sequence,
                timestamp=dt.datetime.now().isoformat(timespec="milliseconds"),
                implementation=implementation,
                temperature=temperature,
                mode=mode,
                iteration=iteration,
                construct_seconds=construct,
                solve_seconds=solve,
                total_seconds=total,
                fixture=args.fixture,
                note=note,
            ))
            sequence += 1

      if (iteration + 1) % max(1, args.runs // 10) == 0:
        print(f"completed {iteration + 1}/{args.runs} iterations")
  finally:
    worker.close()

  write_csv(rows, csv_path)
  write_summary(rows, summary_path)
  if not args.no_plot:
    write_plots(rows, output_prefix)

  print(f"wrote {csv_path}")
  print(f"wrote {summary_path}")
  if not args.no_plot:
    print(f"wrote {output_prefix.with_name(output_prefix.name + '_cold.png')}")
    print(f"wrote {output_prefix.with_name(output_prefix.name + '_warm.png')}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
