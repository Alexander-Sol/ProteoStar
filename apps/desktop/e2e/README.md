# MsViewer E2E / GUI testing (Playwright)

This lets an agent (or a human) drive and inspect the MsViewer GUI with a
Playwright-style API. It's built on
[`tauri-plugin-playwright`](https://github.com/srsholmes/tauri-playwright) —
Tauri apps use the system webview (WebView2 on Windows), which plain Playwright
can't reach, so the plugin embeds a small control server in the app and the test
runner sends it commands.

## Two modes

| Mode      | Folder        | Needs a Rust build? | What it drives                          |
| --------- | ------------- | ------------------- | --------------------------------------- |
| `browser` | `e2e/browser` | No                  | Headless Chromium, Tauri IPC **mocked** |
| `tauri`   | `e2e/tauri`   | Yes (`e2e-testing`) | The **real** MsViewer webview           |

Use **browser** mode for fast, pure-frontend checks (layout, component logic
with mocked backend responses). Use **tauri** mode to test the app for real —
opening datasets, plotting, viewer windows.

> **Platform note (Windows):** the plugin's control server listens on TCP
> `127.0.0.1:6274`. The crate's own `tauri` fixture only speaks Unix sockets, so
> `tauri.fixtures.ts` here is a thin TCP replacement built from the crate's
> exported `PluginClient` + `TauriPage` + `tauriExpect`. Override the port with
> `PW_TAURI_PORT` if 6274 is taken.

## One-time setup

```bash
bun install                       # pulls @playwright/test + @srsholmes/tauri-playwright
bun x playwright install chromium # browser mode only
```

## Running — browser mode

No app needed; the config boots vite automatically.

```bash
cd apps/desktop
bun run e2e:browser
```

## Running — tauri mode (real app)

Two terminals. First build+launch the app **with the testing feature** (first
build is slow; it compiles the plugin):

```bash
# terminal 1 — from apps/desktop
bun run e2e:app          # = tauri dev --features e2e-testing
```

Wait for the window to appear and for this line in its log:

```
tauri-plugin-playwright: listening on tcp://127.0.0.1:6274
```

Then run the tauri specs:

```bash
# terminal 2 — from apps/desktop
bun run e2e:tauri
```

The fixture polls the port, so ordering between the two is forgiving. `bun run
e2e` runs both projects (browser will start vite; tauri needs the app already
up).

## Writing tests

- Browser specs import from `../browser.fixtures` — `tauriPage` is a real
  Playwright `Page` (full Playwright API). Add backend stubs via `ipcMocks` in
  `browser.fixtures.ts`.
- Tauri specs import from `../tauri.fixtures` — `tauriPage` is a `TauriPage` with
  a Playwright-like API (`click`, `fill`, `locator`, `getByRole`,
  `waitForSelector`, `listWindows`, `screenshot`, assertions via `expect`).

The `explore:` tests in each folder show a pattern for discovering the UI (dump
visible text, list open windows) before asserting.

### Real-data example (`tauri/open-and-spectrum.spec.ts`)

Opens `SmallCalibratible_Yeast.mzML`, clicks the TIC apex, and asserts the loaded
MS1 spectrum's base (tallest) peak is **m/z 356.1919** — a full open → click →
read-result round trip against the real backend. It relies on:

- **`e2e/data/`** — git-tracked test fixtures (`SmallCalibratible_Yeast.mzML` +
  its `smalldb.fasta`). `dataPath(name)` resolves them to absolute paths.
- **`e2e/tauri.helpers.ts`** — the reusable moves: `openDatasetByPath`,
  `waitForTic`, `readTicApex`, `clickTicAtRt`, `waitForSpectrum`,
  `readSpectrumBasePeak`. They drive the app's real handlers via
  `TauriPage.evaluate` (find button, emit Plotly's `plotly_click`, read the
  graph div's data) — no reliance on `data-testid`.

### The `window.__msviewerOpenPath` test hook

The native "Open file…" dialog can't be automated (it's an OS window, not DOM).
So `App.tsx` exposes `window.__msviewerOpenPath(path)` — the same open flow the
button uses, minus the dialog. It's installed only when
`window.__PW_ACTIVE__` is set (i.e. the plugin's control server is running), so
it's inert in normal/production builds. `openDatasetByPath` calls it.

## How it's wired into the app

- `src-tauri/Cargo.toml` — `e2e-testing` feature adds the optional
  `tauri-plugin-playwright` dep and `tauri/dynamic-acl`.
- `src-tauri/src/main.rs` — under `#[cfg(feature = "e2e-testing")]`, mounts the
  plugin (TCP 6274) and registers the capability at runtime.
- `src-tauri/e2e-capabilities/playwright.json` — grants `playwright:default`
  across all windows. It lives **outside** `capabilities/` on purpose, so normal
  (non-`e2e-testing`) builds never reference a permission that isn't compiled in.

Nothing here touches release builds: without `--features e2e-testing` the plugin
isn't linked and no server starts.

## macOS / Linux

The same setup works: the plugin prefers its Unix socket there. You can point the
crate's stock `createTauriTest({ mode: 'tauri', mcpSocket: '/tmp/tauri-playwright.sock' })`
at it, or extend `tauri.fixtures.ts` to use the socket path instead of TCP.
