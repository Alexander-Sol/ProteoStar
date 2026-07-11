# ProteoStar — merge plan (flashlfq-rust ⨝ MsBrowser)

Goal: create a **new repo `ProteoStar`** that is a clean copy of **flashlfq-rust** (the host)
with **MsBrowser** merged on top, so that one Cargo workspace holds the engine and both of its
frontends — the **Python bindings** (`flashlfq-py`) and the **desktop app** (`msviewer-tauri` +
the TS/Vite UI). The vendored `crates/flashlfq-core` copy is deleted; the desktop shell path-deps
the real workspace core. Written 2026-07-10.

## Ground truth this plan is built on

- **Host repo:** `F:\flashlfq-rust` (github `Alexander-Sol/flashlfq-rust`). Cargo workspace at
  `rust/` with members `flashlfq-core` + `flashlfq-py`. Currently checked out on branch
  **`topdown-feature-detection`** — this is the **newest core** and the reference for
  `trace_kernel.rs` (per decision). Base ProteoStar on THIS branch, not `main`.
- **Guest repo:** `F:\MsBrowser` (github `Alexander-Sol/MsBrowser`). Bun/TS monorepo
  (`apps/desktop`, `apps/web`, `packages/*`) + vendored `crates/flashlfq-core` +
  `apps/desktop/src-tauri` (standalone Tauri crate). `.git` ≈ 145 MB (fine to carry).
- **The bloat problem:** flashlfq-rust `.git` ≈ **4.7 GB**; an untracked **`analysis/` ≈ 11.2 GB**
  of raw data sits in the working tree; there are ~10 branches (`analysis`, `hill-building-detector`,
  `neighbor-aware-refine`, `obo-envelope-probe`, `overnight-todo-sweep`, worktrees…). Tracked size on
  the target branch is only **215 MB**, so most bloat is history/other-branch raw data. **Do not copy
  the folder and do not push full history** — carry only the target branch and strip large blobs.

## Decisions (locked by the user 2026-07-10)

1. **Keep the engine's commit history** (single-branch clone + `git filter-repo`, not fresh-init).
2. **Leave the `analysis/` folder behind entirely** — it's stale, and new analysis workflows will be
   built around the viewer. The original `F:\flashlfq-rust` stays untouched as the permanent archive of
   that raw data, so nothing is lost; ProteoStar simply never carries it (it also can't — GitHub won't
   host multi-GB blobs). Strip any `analysis/` blobs from ProteoStar's history in Phase 1.
3. **Retire `apps/web`** (the Next.js app) — don't carry it into ProteoStar.
4. **Layout convention (from research):** a **single Cargo workspace at the repo root**, library
   crates under `crates/`, and the **Tauri app crate co-located with its frontend** at
   `apps/desktop/src-tauri`, added to the root workspace as a member. Tauri "mostly doesn't care"
   about monorepo shape — its only rule is `tauri.conf.json` next to the app `Cargo.toml` with
   `frontendDist`/`beforeDevCommand` pointing at the frontend. Keeping `src-tauri` beside the frontend
   means **zero `tauri.conf.json` path surgery**; the root workspace unifies the lockfile + `target/`.
   (Sources: [Tauri project-structure docs](https://v2.tauri.app/start/project-structure/),
   [Tauri monorepo discussion #7368](https://github.com/orgs/tauri-apps/discussions/7368),
   [multi-app discussion #13941](https://github.com/orgs/tauri-apps/discussions/13941).)
- **trace_kernel reconciliation is automatic:** we keep flashlfq-rust's `trace_kernel.rs` by simply
  deleting MsBrowser's vendored core and pointing the Tauri crate at `rust/flashlfq-core`. No manual
  merge of the drifted file is required.
- **Version compatibility (already holds):** flashlfq-rust pins `arrow 54 ↔ pyo3 0.23`; the Tauri
  shell is also on `arrow 54`. One workspace can build the py `cdylib` and the Tauri `bin` together.
  Thermo pin `thermorawfilereader = dotnetrawfilereader-sys = 0.7.0` is load-bearing (0.7.1 breaks
  mzdata 0.65) and must be re-asserted in the unified lockfile.

## Prerequisites

- `git-filter-repo` installed (`pip install git-filter-repo`) — used to strip raw-data blobs from
  history. (If you'd rather not preserve history at all, see the "fresh-history" fallback in Phase 1.)
- Rust toolchain + .NET 8 (thermo reader), `bun`, and `maturin` available.
- Pick a target directory for the new repo, e.g. `F:\ProteoStar`.
- **Do all of this against copies. Leave `F:\flashlfq-rust` and `F:\MsBrowser` untouched** until the
  new repo builds green, so you always have a rollback.

## Phase 0 — Pre-flight in the two source repos (non-destructive)

1. In `F:\flashlfq-rust`, confirm the target branch and that its tree is clean/committed:
   ```
   git -C F:\flashlfq-rust status
   git -C F:\flashlfq-rust rev-parse --abbrev-ref HEAD      # expect topdown-feature-detection
   ```
   Commit or stash any wanted working-tree changes first — the clone in Phase 1 only takes commits.
2. In `F:\MsBrowser`, the working tree has **many uncommitted M/?? changes** (the whole MsViewer +
   feature-finder integration). **Commit them on a branch first** — the subtree merge in Phase 3 only
   brings over committed history:
   ```
   git -C F:\MsBrowser switch -c pre-proteostar-snapshot
   git -C F:\MsBrowser add -A
   git -C F:\MsBrowser commit -m "Snapshot MsViewer + feature-finder work before ProteoStar merge"
   ```

## Phase 1 — Create the lean ProteoStar base from flashlfq-rust

**1a. Single-branch clone (drops other branches' objects, incl. the `analysis` raw data):**
```
git clone --single-branch --branch topdown-feature-detection F:\flashlfq-rust F:\ProteoStar
cd F:\ProteoStar
git branch -m main                      # target branch becomes the new main
```
A single-branch clone repacks only objects reachable from that branch, so blobs that live solely on
`analysis`/other branches are excluded automatically.

**1b. Measure the result:**
```
# size of .git after the clone
(Get-ChildItem F:\ProteoStar\.git -Recurse -File | Measure-Object Length -Sum).Sum/1MB
```
- If it's small (a few hundred MB): good, skip 1c.
- If it's still large: raw data was committed on the target branch's own ancestry — run 1c.

**1c. Strip residual large blobs from history — engine history is kept, analysis is dropped
(decision 1 + 2).** Since `analysis/` is being retired, strip it and any large raw-data blobs from
ProteoStar's history regardless of 1b, so the clone can push to GitHub:
```
git -C F:\ProteoStar filter-repo --analyze          # writes a size report under .git/filter-repo/
# drop the analysis dir + raw-data blobs from ALL of ProteoStar's history:
git -C F:\ProteoStar filter-repo --invert-paths ^
  --path analysis/ ^
  --path-glob "*.raw" --path-glob "*.mzML" --path-glob "*.mzMLb" ^
  --path-glob "*.parquet" --path-glob "*.d"
```
Add any other big paths the `--analyze` report flags. This preserves the engine's *commit history*
(messages, code diffs) while discarding the heavy blobs — exactly the "keep history, drop analysis"
outcome. `filter-repo` strips the origin remote by design; we set a new one in Phase 5. Re-measure
`.git` afterward and confirm it's push-sized (well under 1 GB).

> NOTE: verify the parity goldens under `crates/flashlfq-core/parity/` (the engine's regression
> reference) are NOT caught by the `*.parquet`/`*.raw` globs above — check the `--analyze` report and
> exclude those paths from the strip if needed, or they'll vanish from history too.

**1d. Prune known cruft** in the new tree: `F/`, `rust/F/` (empty artifact dirs), `.venv/`,
`ralph-loop.ps1` (drop unless still used), and the leftover `analysis/` working-tree data (untracked;
delete from the ProteoStar copy — it stays safe in `F:\flashlfq-rust`). Commit the pruning.

## Phase 2 — Decide layout & lockfile policy (before merging MsBrowser)

**Target layout (decision 4 — root Cargo workspace, `crates/`, src-tauri co-located):**
```
F:\ProteoStar\
  Cargo.toml            # [workspace] at REPO ROOT: members = crates/*, apps/desktop/src-tauri
  Cargo.lock            # ONE lockfile — COMMITTED (see below). target/ also lands at repo root.
  crates/
    flashlfq-core/      # <- moved from rust/flashlfq-core (engine + parity goldens + tests + examples)
    flashlfq-py/        # <- moved from rust/flashlfq-py  (PyO3 bindings → maturin)   frontend #1
  apps/
    desktop/            # <- MsBrowser apps/desktop (TS/Vite UI)
      src-tauri/        # <- stays here, beside its frontend; a member of the root workspace  frontend #2
  packages/             # <- MsBrowser packages/* (ui, plot-adapter, viewer-state, imsp-core)
  package.json          # <- MsBrowser bun-workspaces root (name → "proteostar")
  tsconfig.base.json, tsconfig.packages.json, vitest.config.ts   # <- MsBrowser root TS config
  agent_info/, README.md, CLAUDE.md                              # merged docs
```
Two moves vs. the old flashlfq-rust layout: the workspace root goes from `rust/` up to the repo root,
and `rust/flashlfq-{core,py}` become `crates/flashlfq-{core,py}`. The **Tauri crate does NOT move** —
it stays at `apps/desktop/src-tauri` (beside `tauri.conf.json` and the frontend, exactly where Tauri
expects it) and is pulled into the root workspace by path. No `tauri.conf.json` path surgery.
`rust/` is deleted once its two crates are relocated.

**Root `Cargo.toml` workspace:**
```toml
[workspace]
resolver = "2"
members = ["crates/flashlfq-core", "crates/flashlfq-py", "apps/desktop/src-tauri"]
default-members = ["crates/flashlfq-core", "apps/desktop/src-tauri"]   # see Phase 4.4 re: py/libpython

[workspace.package]      # carried over from rust/Cargo.toml
version = "0.1.0"
edition = "2021"
license = "MIT"
authors = ["mzLib contributors"]
repository = "https://github.com/Alexander-Sol/ProteoStar"

[workspace.dependencies]
pyo3 = { version = "0.23", default-features = false }
numpy = "0.23"
arrow = "54"
rayon = "1.10"
serde = { version = "1", features = ["derive"] }
```

**Lockfile policy change:** flashlfq-rust's `.gitignore` currently ignores `Cargo.lock` (library
convention). ProteoStar ships **binaries** (the Tauri app + the py extension), so **commit the root
workspace `Cargo.lock`** — this is where the thermo=0.7.0 pin is preserved. Remove the `Cargo.lock`
line from the merged `.gitignore`, and delete the two now-stale nested lockfiles
(`apps/desktop/src-tauri/Cargo.lock`, the vendored core's lock).

## Phase 3 — Merge MsBrowser history in (subtree, history preserved)

From `F:\ProteoStar` (on `main`), pull MsBrowser's committed snapshot in under a temporary prefix,
then move files into the target layout so history follows the moves:
```
git remote add msbrowser F:\MsBrowser
git fetch msbrowser
# bring the whole MsBrowser tree in under _msbrowser/ preserving history:
git merge -s ours --no-commit --allow-unrelated-histories msbrowser/pre-proteostar-snapshot
git read-tree --prefix=_msbrowser/ -u msbrowser/pre-proteostar-snapshot
git commit -m "Merge MsBrowser history under _msbrowser/ (temporary staging)"
```
Before this, relocate flashlfq-rust's own crates to the new layout (also history-preserving):
```
git mv rust/flashlfq-core crates/flashlfq-core
git mv rust/flashlfq-py   crates/flashlfq-py
# move the root workspace Cargo.toml up and rewrite members/paths (Phase 2), then:
git rm -r rust            # now empty except the old workspace files you've relocated
```
Then bring MsBrowser's tree into the layout:
```
git mv _msbrowser/apps/desktop apps/desktop     # frontend + its src-tauri (stays put)
git mv _msbrowser/packages packages
git mv _msbrowser/package.json _msbrowser/tsconfig.base.json _msbrowser/tsconfig.packages.json _msbrowser/vitest.config.ts .
git mv _msbrowser/bun.lock .
# merge docs rather than clobber: move MsBrowser's *.md / agent_info into agent_info/ or docs/
# DELETE the vendored core — the real one at crates/flashlfq-core wins (this is what keeps
# flashlfq-rust's trace_kernel.rs as canonical, per the decision):
git rm -r _msbrowser/crates/flashlfq-core
# RETIRE apps/web (decision 3) and drop other cruft (.DS_Store, screenshots zip, etc.):
git rm -r _msbrowser/apps/web
rmdir _msbrowser      # should be empty now; commit
```
Result: both projects' histories are in one repo, files carry blame through the moves.

## Phase 4 — Unify the Cargo workspace

1. **Root workspace `Cargo.toml`** as written in Phase 2 — members `crates/flashlfq-core`,
   `crates/flashlfq-py`, `apps/desktop/src-tauri`. The Tauri crate joins by path; **it does not move**.
2. **Repoint the Tauri crate's core dep** in `apps/desktop/src-tauri/Cargo.toml`:
   `flashlfq-core = { path = "../../../crates/flashlfq-core" }`  (three up from `apps/desktop/src-tauri`
   to repo root, then into `crates/`). Verify the `..` count against the final tree.
3. **Hoist shared deps** to `[workspace.dependencies]` where they overlap (`arrow`, `serde`, `rayon`)
   and reference them with `arrow = { workspace = true }` in each crate, so core/py/tauri stay in
   lockstep. Update the `flashlfq-core` `Cargo.toml` header comment — it's back in a workspace, so the
   "de-workspaced pins" note no longer applies.
4. **Keep `flashlfq-py` out of `default-members` if a bare `cargo build` fails to link libpython.**
   The Phase 2 `Cargo.toml` already sets `default-members = ["crates/flashlfq-core",
   "apps/desktop/src-tauri"]` so `cargo build`/`cargo test` at the root build core+tauri and the wheel
   is built explicitly via `maturin`. Confirm this is actually needed on this machine (it usually is
   for `extension-module` crates) and keep or drop accordingly.
5. **Re-assert the thermo pin** in the unified root `Cargo.lock`:
   ```
   cargo update -p thermorawfilereader --precise 0.7.0
   cargo update -p dotnetrawfilereader-sys --precise 0.7.0
   ```
6. **No `tauri.conf.json` path changes** — because `src-tauri` stayed beside its frontend, the
   `beforeDevCommand`/`frontendDist`/`devUrl` in `apps/desktop/src-tauri/tauri.conf.json` are already
   correct. The only new effect is that `target/` now resolves to the **workspace root** (shared by
   all crates) instead of per-crate; confirm nothing hard-codes `apps/desktop/src-tauri/target`.

## Phase 5 — New identity & remote

1. Merge `.gitignore`s (union of both; **remove** the `Cargo.lock` ignore per Phase 2). Merge
   `README.md` / `CLAUDE.md` into ProteoStar-branded versions. Update `package.json` `name` from
   `msbrowser` → `proteostar`.
2. Create the empty GitHub repo `ProteoStar` (via `gh repo create Alexander-Sol/ProteoStar --private`
   or the web UI). Then:
   ```
   git -C F:\ProteoStar remote remove origin 2>$null
   git -C F:\ProteoStar remote remove msbrowser
   git -C F:\ProteoStar remote add origin https://github.com/Alexander-Sol/ProteoStar.git
   git -C F:\ProteoStar push -u origin main
   ```
   Confirm the push size is sane (the whole point of Phase 1 — if it tries to push multiple GB, stop
   and revisit 1c). GitHub warns >1 GB and rejects individual files >100 MB.

## Phase 6 — Verify both frontends build off the one core

- Engine unchanged: `cargo test -p flashlfq-core` (runs the restored parity goldens).
- Desktop backend: `cargo build -p msviewer-tauri`.
- Python frontend: `cd crates/flashlfq-py && maturin develop` (or `maturin build`), then run a smoke
  test (`python quant_smoke_test.py` / `xic_smoke_test.py`).
- Full desktop app: `cd apps/desktop && bun install && bun run tauri dev --release` and drive it on the
  Jurkat raw (`D:\JurkatTopdown\02-18-20_jurkat_td_rep1_fract6.raw`) — open dataset, run feature
  finding, confirm the overlay still works end-to-end.
- TS unit tests: `bun run test` at root (vitest).

Only after all six are green: delete/park the old `F:\MsBrowser` and `F:\flashlfq-rust` working copies
(keep them until you're confident — they're the rollback).

## Phase 7 — Follow-ups (not blocking the merge)

- **CI:** neither source repo has `.github/workflows`. Add one ProteoStar pipeline: cargo
  build+test (core), maturin build (py wheel), bun build+test + tauri build (desktop).
- **Update the handoff docs** (`MsViewer_Architecture.md` §4/§12.7, `FeatureFinder-Integration.md`
  "re-vendor" limitation) to say the core is now a workspace sibling, not a vendored snapshot — the
  whole re-sync problem is gone.
- **Prune stale flashlfq-rust branches** you don't carry into ProteoStar (they're only in the old
  repo now).

## Decisions — all locked (2026-07-10)

1. ✅ **Keep engine history** — single-branch clone + `git filter-repo` (strips analysis/raw blobs,
   preserves commit history).
2. ✅ **Drop `analysis/`** — stays only in the untouched `F:\flashlfq-rust`; never enters ProteoStar.
3. ✅ **Retire `apps/web`** — `git rm -r` during the merge (Phase 3).
4. ✅ **Layout** — root Cargo workspace, `crates/flashlfq-{core,py}`, `apps/desktop/src-tauri`
   co-located and joined by path (Phase 2/4). No `tauri.conf.json` surgery.

Nothing blocks execution — ready to run Phase 0 on the user's go-ahead.
