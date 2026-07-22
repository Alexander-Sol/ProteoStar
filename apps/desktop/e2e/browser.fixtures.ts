// Browser-mode fixtures: headless Chromium with Tauri IPC mocked in JS.
// No Rust build, no running app — the vite dev server (started by the config's
// webServer) serves the frontend and `createTauriTest` shims window.__TAURI__.
//
// `tauriPage` here is a real Playwright Page, so the full Playwright API is
// available. Add ipcMocks below as tests start exercising commands that hit the
// backend (open_dataset, get_spectrum, …) — unmapped commands resolve undefined.
import { createTauriTest } from "@srsholmes/tauri-playwright";

export const { test, expect } = createTauriTest({
  devUrl: "http://localhost:1420",
  ipcMocks: {
    // MsViewer calls the backend only after a dataset is opened, so an empty
    // map is enough for render/smoke tests. Extend per-test as needed, e.g.:
    // get_metadata: () => ({ scanCount: 0, ... }),
  },
});
