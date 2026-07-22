import { test, expect } from "../browser.fixtures";

// Browser-mode smoke test: proves the harness renders MsViewer's frontend in
// headless Chromium with Tauri IPC mocked. `tauriPage` is a real Playwright Page.
test("app shell renders", async ({ tauriPage }) => {
  await expect(tauriPage).toHaveTitle(/MsViewer/i);
  // React mounts into #root — a non-empty root means the app booted.
  await expect(tauriPage.locator("#root")).not.toBeEmpty();
});

test("explore: dump the initial visible text", async ({ tauriPage }) => {
  // A pattern agents can copy to discover what's on screen before asserting.
  const text = await tauriPage.locator("body").innerText();
  console.log("[MsViewer visible text]\n" + text.slice(0, 2000));
  expect(text.length).toBeGreaterThan(0);
});

// Task 1 (no-scroll layout), fast path: the shell must not overflow the viewport vertically.
// Browser mode loads no spectrum, so this guards the empty-shell layout without a Rust build;
// the real spectrum-loaded case is covered by the tauri spec.
test("the shell does not scroll vertically", async ({ tauriPage }) => {
  const scrollHeight = await tauriPage.evaluate<number>("document.documentElement.scrollHeight");
  const innerHeight = await tauriPage.evaluate<number>("window.innerHeight");
  console.log(`[layout] scrollHeight=${scrollHeight} innerHeight=${innerHeight}`);
  expect(scrollHeight).toBeLessThanOrEqual(innerHeight + 2);
});
