"""Run the existing Marlpde scenario with LSODA instrumentation.

Use the model's Pipenv environment; this analysis environment lacks py-pde.
No model files are modified. Output is diagnostic, not manuscript input.
"""

import argparse
from dataclasses import asdict
import faulthandler
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.integrate import LSODA, solve_ivp


def diagnostic_solver(directory, tstar, interval=10):
    """Construct an instrumented LSODA class for solve_ivp."""
    class DiagnosticLSODA(LSODA):
        def __init__(self, fun, *args, **kwargs):
            self.started = self.reported = time.monotonic()
            self.accepted_steps = self.calls = 0
            self.last_h = 0.0

            def checked(t, y):
                self.calls += 1
                self.check("state", t, y)
                if time.monotonic() - self.reported >= interval:
                    self.report("rhs", t)
                result = fun(t, y)
                self.check("rhs", t, result)
                return result

            super().__init__(checked, *args, **kwargs)
            self.check("initial", self.t, self.y)
            self.checkpoint()
            self.report("initial", self.t)

        def check(self, label, t, values):
            if not np.isfinite(t) or not np.isfinite(values).all():
                np.savez(directory / "invalid_evaluation.npz", t=t,
                         values=values, kind=label)
                raise FloatingPointError(
                    f"Non-finite {label} at trial t={t:.17g}; "
                    "see invalid_evaluation.npz and checkpoint.npz"
                )

        def checkpoint(self):
            temporary = directory / "checkpoint.tmp.npz"
            np.savez(temporary, t=self.t, y=self.y, Tstar=tstar,
                     accepted_steps=self.accepted_steps)
            temporary.replace(directory / "checkpoint.npz")

        def report(self, kind, trial_t):
            self.reported = time.monotonic()
            record = dict(kind=kind, wall_s=self.reported-self.started,
                          accepted_t=float(self.t), accepted_kyr=self.t*tstar/1000,
                          trial_t=float(trial_t), last_h=self.last_h,
                          steps=self.accepted_steps, rhs_calls=self.calls,
                          njev=int(self.njev), nlu=int(self.nlu),
                          minimum=float(np.min(self.y)), maximum=float(np.max(self.y)))
            with (directory / "progress.jsonl").open("a") as stream:
                stream.write(json.dumps(record) + "\n")
            print(json.dumps(record), flush=True)
            self.checkpoint()

        def step(self):
            before = self.t
            try:
                message = super().step()
                self.check("accepted", self.t, self.y)
                if self.status != "failed":
                    self.accepted_steps += 1
                    self.last_h = float(self.t-before)
                if self.status != "running" or time.monotonic()-self.reported >= interval:
                    self.report(self.status, self.t)
                return message
            except BaseException:
                # During an interrupted step, the on-disk checkpoint from the
                # previous report remains the known accepted state.
                raise

    return DiagnosticLSODA


class DiagnosticFinished(Exception):
    pass


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-directory", type=Path, required=True,
                        help="directory containing Evolve_scenario.py")
    parser.add_argument("--output-directory", type=Path, required=True,
                        help="new directory outside the analysis repository")
    parser.add_argument("--end-time", type=float,
                        help="optional final dimensionless time for a short smoke test")
    args = parser.parse_args()
    sys.path.insert(0, str(args.model_directory.resolve()))
    import Evolve_scenario as model

    pde = asdict(model.Map_Scenario())
    solver = asdict(model.Solver())
    tracker = asdict(model.Tracker())
    if solver["method"] != "LSODA":
        parser.error("Set the model's Solver.method to LSODA first")
    if solver.get("lband") is not None or solver.get("uband") is not None:
        parser.error("Set lband and uband to None first")
    if args.end_time is not None:
        start, end = solver["t_span"]
        if not start < args.end_time <= end:
            parser.error("--end-time must be inside the configured integration interval")
        solver["t_span"] = (start, args.end_time)
        tracker["t_eval"] = np.linspace(start, args.end_time, 2)
    directory = args.output_directory.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    (directory / "parameters.json").write_text(json.dumps(
        dict(model_directory=str(args.model_directory.resolve()),
             model=pde, solver=solver, tracker=tracker), indent=2,
        default=lambda value: value.tolist() if isinstance(value, np.ndarray) else str(value)))

    # Keep model construction and event handling, replace only its solver call.
    def run(fun, **options):
        options["method"] = diagnostic_solver(directory, float(pde["Tstar"]))
        result = solve_ivp(fun, **options)
        np.savez(directory / "sampled_solution.npz", t=result.t, y=result.y,
                 Tstar=pde["Tstar"], success=result.success, message=result.message)
        print(result.message, flush=True)
        if not result.success:
            raise RuntimeError(result.message)
        # Bypass the model's normal storage and plotting, including its paths.
        raise DiagnosticFinished

    model.solve_ivp = run
    # Disable the old trial-time progress bar, which is not accepted progress.
    original_tqdm = model.tqdm
    model.tqdm = lambda *a, **kw: original_tqdm(*a, **kw, disable=True)
    print(f"Diagnostics: {directory}", flush=True)
    with (directory / "stacks.log").open("w") as stacks:
        faulthandler.dump_traceback_later(60, repeat=True, file=stacks)
        try:
            model.integrate_equations(solver, tracker, pde)
        except DiagnosticFinished:
            pass
        finally:
            faulthandler.cancel_dump_traceback_later()


if __name__ == "__main__":
    main()
