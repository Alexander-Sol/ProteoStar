import { defineConfig, devices } from "@playwright/test";

// Two ways to drive MsViewer, split by folder so each has its own transport:
//
//   browser  → e2e/browser/*.spec.ts — headless Chromium with Tauri IPC mocked.
//              Fast, needs no Rust build. `webServer` boots vite for it.
//   tauri    → e2e/tauri/*.spec.ts   — the REAL webview via the plugin's control
//              server (TCP 127.0.0.1:6274 on Windows). Requires the app running
//              with `bun run e2e:app` (tauri dev --features e2e-testing).
//
// See e2e/README.md for the full workflow.
export default defineConfig({
  fullyParallel: false,
  forbidOnly: !!process.env.CI,
  retries: 0,
  reporter: [["list"], ["html", { open: "never", outputFolder: ".report" }]],
  timeout: 30_000,
  outputDir: ".artifacts",

  projects: [
    {
      name: "browser",
      testDir: "./browser",
      use: {
        ...devices["Desktop Chrome"],
        trace: "on",
        screenshot: "only-on-failure",
      },
    },
    {
      name: "tauri",
      testDir: "./tauri",
      // Traces/screenshots from Playwright's own browser are meaningless here —
      // the fixture talks to the native webview over a socket. The tauri fixture
      // captures native screenshots on failure instead.
      use: { trace: "off", screenshot: "off" },
    },
  ],

  // Only the browser project needs the vite dev server. The tauri project's app
  // is launched separately (bun run e2e:app), so its fixture just connects.
  webServer: process.env.PW_NO_WEBSERVER
    ? undefined
    : {
        command: "bun x vite --port 1420",
        port: 1420,
        cwd: "..",
        reuseExistingServer: true,
        timeout: 60_000,
      },
});
