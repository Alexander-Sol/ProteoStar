# Top-Down Simulator — persistent live session

**As of:** 2026-10-20 · **Branch:** `sim-workbench-design` · precedes Phase 3 (the workbench UI) of
`TopDown-Simulator-Rust-Port-Plan.md`.

## Why

The live path rebuilds a `Simulator` from the full spec on every edit. Measured with
`examples/simulate_timing` (1000 proteoforms × 3000 scans, charges 3–60, default noise, release):

| Step | Time |
|---|---|
| `Simulator::new` (paid on every `simulate_scan` call today) | 73 ms |
| `simulate_scan` on that instance | 29 ms |
| First `simulate_scan` with **noise off** | ~28 s: it centroids the whole run to find the brightest centroid |

With noise on, ~70 % of each edit is rebuilding envelopes the edit did not touch. With noise off
the live path is unusable at run size. The fix is one long-lived `Simulator` per workbench that
takes each new spec as an **update** and recomputes only what the update invalidates.

## Design

### 1. `Simulator::update(&mut self, spec)` — diff, invalidate, keep the rest

The simulator holds the spec it currently simulates and diffs each new one against it. Rule: after
`update(s)`, every output is **bit-identical** to a fresh `Simulator::new(s)`. Caching is never
allowed to change results; that is the test.

What each edit invalidates:

| Edit | Recomputed | Kept |
|---|---|---|
| Abundance, RT μ/σ/τ | that proteoform's model (envelopes don't depend on it) | all envelopes |
| Charge μ_z/σ_z | that proteoform's envelopes, from its cached isotope kernel | other envelopes, kernel |
| Mass | that proteoform's kernel + envelopes | other envelopes |
| Add / remove proteoform | its envelopes (spliced in/out) | other envelopes |
| Width `k`, charge range | every envelope, from the cached kernels | kernels |
| Scan times | scan-indexed caches | envelopes |
| Noise, reduction | nothing in the signal model | envelopes, scan maxima, viewed-scan signal (if its prune floor is unchanged) |

Proteoform diff = longest common prefix + suffix of the combined `proteoforms ++ background` list.
The changed middle is replaced. When the lengths match it is compared element by element, so a
single-field edit rebuilds at most one proteoform. Loading an unrelated spec makes the whole list
the middle, which costs the same as `new`.

Structural changes in `ProfileCentroider`:
- Envelopes are stored per proteoform (`Vec<Vec<Envelope>>`), so add/remove is a splice. Iteration
  order equals the old flat list's, so outputs stay bit-identical to the goldens' pipeline.
- Each proteoform's `IsotopeEnvelope` kernel is cached, so `k` / charge-range edits skip averagine.
- Every fallible step (envelope building) runs before anything is committed, so a failed update
  leaves the simulator unchanged.

### 2. The brightest centroid, per scan, cached

Without noise, the reduction floor is `relative × the brightest centroid in the run`. It depends
only on the signal model (proteoforms, `k`, charge range, scan times), not on noise or reduction.

- **Store one maximum per scan** (`Vec<Option<ScanMax>>`, with `None` meaning stale), not one
  number. Toggling noise, or editing noise/reduction, never recomputes it.
- **Compute a scan's maximum cheaply.** Centroid the scan with its own tallest Gaussian `T` as the
  floor bound (`centroid_above(t, T)`). The brightest centroid is ≥ `T`, so it is always
  reported. Everything below 1e-6·`T` is pruned, so only the neighbourhood of the giants is walked.
- **Invalidate locally.** An edit changes a set of proteoforms. Scan `s` is stale only if one of
  them, before or after the edit, contributes a Gaussian at or above `1e-6·T_s`, the scan's prune
  threshold. Below it the scan's component set and its `T` are provably identical, so the cached
  value is exactly what a fresh build would compute. An EMG is localised in RT, so editing one
  proteoform marks tens of scans stale, not thousands. `k`, charge-range and scan-time edits mark
  every scan stale.
- **Fallback if needed.** If refreshing every stale scan ever proves too slow for the live view,
  preview can use the cached maxima plus the tallest-Gaussian bound for stale scans, and stay exact
  on Accept. Measure first.

Consequence: the noise-free floor is known before centroiding, so scans are pruned against the
floor itself rather than against a looser bound. `new` no longer scans the run for the tallest
Gaussian.

### 3. Viewed-scan signal cache

The pre-noise, pre-reduction centroids of the last scan simulated are kept with the prune floor
they were computed at. A noise-only edit re-thresholds and re-injects them without centroiding.
The prune floor is **quantised down to a power of two** before use, so small noise edits keep the
same prune floor. That keeps the cache valid while outputs stay a pure function of the spec. Any
signal edit clears this cache.

### 4. Tauri session commands (`src-tauri`)

`AppState.sim_sessions: HashMap<u64, Arc<Mutex<Simulator>>>`.

- `open_sim_session(spec) → session`
- `sim_session_scan(session, spec, scanIndex, mzMin?, mzMax?, maxPeaks?) → spectrum`: update, then
  simulate one scan. One round trip per edit.
- `commit_sim_session(session, spec) → OpenResult`: update, then `simulate_run` on the same
  instance, then index and register it as a dataset (shared with `open_simulated_dataset`).
- `close_sim_session(session)`

The frontend always sends the whole `SimSpec` (panel state = spec). The backend does the diffing.
The stateless `simulate_scan` / `open_simulated_dataset` stay for one-off use.

### 5. Frontend: latest-wins, not debounce

A `SimSession` wrapper in `tauri-dataset-provider.ts` keeps at most one scan request in flight.
Edits during it replace a single pending request, and superseded requests resolve to `null`.
Slider drags then skip intermediate positions with no tuning constant.

## Not now: incremental Accept

With noise on, an edit to one proteoform only changes the scans in its RT window, so Accept could
re-simulate those and patch the committed dataset. It costs the run's pre-noise signal in memory
(hundreds of MB at run size) plus index patching. Revisit if Accept latency becomes the bottleneck.

## Status (2026-10-20): implemented, uncommitted

Items 1–5 are built: `Simulator::update` / `global_max` / `cached_scan_maxima` (`sim-core/src/lib.rs`),
the per-proteoform `Species` + `ProfileCentroider::update` + `scan_max` (`centroider.rs`), the four
session commands (`src-tauri/src/commands.rs`, `state.rs`), and `SimSession` in
`tauri-dataset-provider.ts`. Tests: `sim-core/tests/session_update.rs` (update ≡ new bit-for-bit over
22 chained edits from both a noise-on and a noise-off start, cache survival, and a failed update
leaving state unchanged). A mutation check showed it catches a stale scan-maximum cache. All 23
sim-core tests are green, `msviewer-tauri` checks, and `tsc` is clean.

Measured (`simulate_timing`, 1000 proteoforms × 3000 scans, release, update + viewed scan):

| Edit | Before | Now |
|---|---|---|
| Proteoform edit, noise on | ~100 ms | 28–35 ms (the viewed scan's own centroiding) |
| Noise edit | ~100 ms | 4–5 ms |
| First view with noise off | ~28 s | 1.8 s |
| Proteoform edit, noise off | ~28 s | 80–136 ms (~89 % of scan maxima kept) |
| Noise turned back off | ~28 s | 11 ms |
| `k` edit, noise off | ~28 s | 2.0 s (every scan maximum recomputed) |
| `simulate_run` (Accept) | 28.5 s | 14.8 s (prunes against the exact floor, not a bound) |

The stale-scan fallback (design §2) was not needed. A local edit invalidates ~11 % of scans: the
1e-6 prune threshold reaches well into the EMG tails, and further where `T` is small at the gradient
ends.

## Verification

- `update` ≡ `new`, bit for bit, for each edit kind (abundance, RT, charge, mass, add, remove,
  `k`, charge range, scan times, noise on/off, reduction), on `simulate_scan` and `simulate_run`.
- Goldens and `fast_path_exactness` stay green.
- `simulate_timing` gains session numbers: an update+scan per edit kind, with noise on and off.
