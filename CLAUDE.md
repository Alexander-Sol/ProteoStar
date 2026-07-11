# ProteoStar — project guidance

Unified workspace: the FlashLFQ Rust engine + untargeted MS1 feature detector, its Python
bindings (`crates/flashlfq-py`), and the MsViewer Tauri desktop app (`apps/desktop`). Core crate:
`crates/flashlfq-core`. Detector runner: `crates/flashlfq-core/examples/detect_features_tsv.rs`.
Design docs live in `agent_info/` (start with `Feature-Detection-Design.md`,
`Detector-Perf-and-Parallelization.md`, and `MsViewer_Architecture.md`).

## Benchmarking

**Whenever you benchmark the detector, follow `agent_info/Benchmarking-Guide.md`.** It is the
authoritative procedure — detector invocation → ground-truth construction → recall scoring — and holds
the full input paths (raw + PSM reference), the per-file RT tolerances and reference-peak counts, and
the current default numbers. Do not hand-roll a different benchmark or reference set; defer to the guide.

**Record every benchmark result in the guide** (its §5 Results table): update it with the numbers, date,
and commit so it always holds the most recent results.

Efficiency rules:
- **Do not re-run the default pipeline on these files just to get a comparison number.** The default
  detect time / feature count / recall are already in the guide — diff against those.
- When measuring a change, run **only the changed (non-default) configuration** and diff it against the
  guide's saved baseline. Re-running the default alongside it wastes minutes on the big files.
- Re-run the default path **only when explicitly asked, or when the default itself changes.** In the
  latter case, update the guide's Results table (numbers, date, commit) so it stays the current reference.
- Detect dominates wall-clock on the big files; use `DETECT_ONLY=1` when only detect cost matters. The
  shipped default detect path is 2-D tiling at `COVERAGE_TARGET=1.0` (any coverage cap < 1.0 forces the
  serial fallback). The baseline recall numbers require `APEX_PREGATE=0` until that flag's default is
  made opt-in (see the guide's caveat).
