import { test, expect } from "../tauri.fixtures";

// Tauri-mode smoke test: drives the REAL MsViewer webview through the plugin's
// control server. Requires the app running: `bun run e2e:app`.
test("real webview reports the app title", async ({ tauriPage }) => {
  const title = await tauriPage.title();
  expect(title).toMatch(/MsViewer/i);
});

test("real webview has mounted the React root", async ({ tauriPage }) => {
  await tauriPage.waitForSelector("#root");
  const html = await tauriPage.innerHTML("#root");
  expect(html.length).toBeGreaterThan(0);
});

test("explore: list windows and dump visible text", async ({ tauriPage }) => {
  // Discover every open webview (MsViewer may spawn viewer windows at runtime).
  const windows = await tauriPage.listWindows();
  console.log("[open windows]", JSON.stringify(windows, null, 2));

  const text = await tauriPage.innerText("body");
  console.log("[MsViewer visible text]\n" + text.slice(0, 2000));
  expect(windows.length).toBeGreaterThan(0);
});
