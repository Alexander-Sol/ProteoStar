# ProteoStar — project guidance

Unified workspace: the FlashLFQ Rust engine + untargeted MS1 feature detector, its Python
bindings (`crates/flashlfq-py`), and the MsViewer Tauri desktop app (`apps/desktop`). Core crate:
`crates/flashlfq-core`. Detector runner: `crates/flashlfq-core/examples/detect_features_tsv.rs`.
Design docs live in `agent_info/` (start with `Feature-Detection-Design.md`,
`Detector-Perf-and-Parallelization.md`, and `MsViewer_Architecture.md`).

## Querying the codebase — use serena first

**Use the serena MCP tools for codebase queries.** Serena provides semantic, symbol-aware search
(`find_symbol`, `find_referencing_symbols`, `get_symbols_overview`, etc.) that is more efficient and
precise than raw text search for navigating this workspace. Reach for serena before falling back to
Grep/Glob/Read when you need to locate definitions, callers, or the shape of a module. Per serena's
own instruction, call `initial_instructions` to load its manual before starting a coding task.

## GUI / E2E testing (Playwright) — `apps/desktop/e2e`

The MsViewer desktop app has a Playwright harness for driving and inspecting the GUI. Full docs
live in **`apps/desktop/e2e/README.md`** — read it before writing or running tests. Built on
[`tauri-plugin-playwright`](https://github.com/srsholmes/tauri-playwright) (Tauri apps use the
system webview, which plain Playwright can't reach).

**Two modes, split by folder:**
- `e2e/browser/*.spec.ts` — headless Chromium, Tauri IPC **mocked**. No Rust build. Fast UI/logic
  checks. `tauriPage` is a real Playwright `Page`. Fixture: `e2e/browser.fixtures.ts`.
- `e2e/tauri/*.spec.ts` — the **real** webview + Rust backend. `tauriPage` is a `TauriPage`
  (Playwright-like API). Fixture: `e2e/tauri.fixtures.ts`. This is the one that exercises real data.

**Running (from `apps/desktop`):**
- Browser: `bun run e2e:browser` (boots vite itself).
- Tauri: `bun run e2e:app` in one terminal (= `tauri dev --features e2e-testing`; wait for
  `listening on tcp://127.0.0.1:6274`), then `$env:PW_NO_WEBSERVER=1; bun run e2e:tauri` in another.

**Windows specifics (this is a Windows box):**
- The plugin's control server listens on **TCP `127.0.0.1:6274`** (its Unix socket is unix-only).
  The crate's stock `tauri` fixture is socket-only, so `tauri.fixtures.ts` is a thin TCP fixture
  built from the crate's exported `PluginClient` + `TauriPage` + `tauriExpect`.
- Set `PW_NO_WEBSERVER=1` for tauri runs so Playwright doesn't fight `tauri dev`'s vite for port 1420.
- If a run fails with "Port 1420 already in use" or 6274 busy, kill leftovers:
  `Get-NetTCPConnection -LocalPort 1420,6274 -State Listen | %{ Stop-Process -Id $_.OwningProcess -Force }`
  and `Get-Process msviewer-tauri | Stop-Process -Force`.

**How it's wired (all behind the `e2e-testing` cargo feature; release builds are untouched):**
- `src-tauri/Cargo.toml` — optional `tauri-plugin-playwright` + `tauri/dynamic-acl`.
- `src-tauri/src/main.rs` — under `#[cfg(feature = "e2e-testing")]`, mounts the plugin (TCP 6274)
  and registers the capability at runtime via `app.add_capability(include_str!(...))`.
- `src-tauri/e2e-capabilities/playwright.json` — grants `playwright:default`; kept **outside**
  `capabilities/` so normal builds don't reference an un-compiled permission.

**Writing tauri-mode tests — key facts and gotchas:**
- `TauriPage` selectors go through plain `document.querySelector` — **no** Playwright `:has-text`,
  and `getByText` (maps to a `[data-pw-text]` attr the socket server doesn't stamp) is unreliable.
  The app has no `data-testid`s yet; drive via `TauriPage.evaluate("<js>")` (find button by text,
  read/emit on DOM) — see `e2e/tauri.helpers.ts`.
- **The native "Open file…" dialog cannot be automated.** `App.tsx` exposes an E2E-only hook
  `window.__msviewerOpenPath(path)` (same open flow, no dialog) that exists **only** when
  `window.__PW_ACTIVE__` is set by the plugin — inert in production. `openDatasetByPath` uses it.
  Do NOT try to intercept `__TAURI_INTERNALS__.invoke` — that property isn't writable, and clicking
  the real "Open file…" button opens the native picker and **freezes the webview**.
- Plots are Plotly; helpers fire the app's real handlers with `gd.emit('plotly_click', …)` and read
  results off the graph div's `.data`/`.layout` (e.g. `readSpectrumBasePeak`). No pixel math.
- Each `TauriPage.evaluate` is capped at **30s** server-side — poll from Node across multiple
  evaluates, don't loop inside one.
- Test data lives in `e2e/data/` (git-tracked fixtures); `dataPath(name)` → absolute path.
- Reusable moves in `e2e/tauri.helpers.ts`: `openDatasetByPath`, `waitForTic`, `readTicApex`,
  `clickTicAtRt`, `waitForSpectrum`, `readSpectrumBasePeak`. Reference test:
  `e2e/tauri/open-and-spectrum.spec.ts` (open mzML → click TIC apex → assert base peak m/z 356.1919).

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
