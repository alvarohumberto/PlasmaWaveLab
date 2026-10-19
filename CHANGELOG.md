# Changelog

## 0.3.0 (2026-10-19)
- New: `critical_mach()` returns the upper Mach limit M_c = 1.5852010652 of the
  soliton branch; `peak_potential` now reports it when M >= M_c.
- New: `ion_acoustic_kdv_coefficients()` and `small_amplitude_soliton()` link
  the Sagdeev and KdV modules through the Washimi–Taniuti small-amplitude limit.
- New: `examples/verification.py` reproduces every figure and number in the
  SoftwareX manuscript (spatial/temporal convergence, Mach sweep, KdV limit,
  two-soliton phase shifts, invariant drift) and writes `summary.json`.
- CLI: new options `--xi-max`, `--length`, `--dt`, `--amplitude`, `--alpha`,
  `--beta`, `--version`; KdV error now reported in the maximum norm.
- Tests: 11 new verification tests, including two-soliton phase shifts (21 in total); CI also runs the
  verification script.
- Manuscript moved to Elsevier `elsarticle` with SoftwareX OSP structure.

## 0.2.0 (2026-10-19)
- GitHub-ready Python project with packaging metadata, CLI entry point, CI matrix, citation and contributor guidance.
- Baseline Sagdeev and KdV solvers, numerical tests and reproducible output generators.
- Model assumptions and limitations documented.

## 0.1.0 (2026-10-19)
- Initial working research prototype.
