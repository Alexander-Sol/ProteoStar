// Tauri-mode fixtures: drive the REAL MsViewer webview through the plugin's
// control server. On Windows the plugin listens on TCP 127.0.0.1:6274 (the Unix
// socket path in the crate is #[cfg(unix)] only).
//
// The crate ships a `createTauriTest({ mode: 'tauri' })` fixture, but it connects
// over a Unix socket path and never wires up TCP — so it can't reach the server
// on Windows. We rebuild the tiny bit that's missing from the crate's own
// exported primitives (PluginClient + TauriPage + tauriExpect): connect over TCP,
// hand the test a TauriPage. Everything downstream (locators, assertions,
// waitForWindow, native screenshots) is the crate's real implementation.
import { test as base } from "@playwright/test";
import {
  PluginClient,
  TauriPage,
  tauriExpect,
} from "@srsholmes/tauri-playwright";
import { createConnection } from "node:net";

const TCP_PORT = Number(process.env.PW_TAURI_PORT ?? 6274);
const HOST = "127.0.0.1";

/** Poll the control-server port until the app is up (it takes a moment to boot). */
async function waitForPort(port: number, timeoutMs = 60_000): Promise<void> {
  const deadline = Date.now() + timeoutMs;
  // eslint-disable-next-line no-constant-condition
  while (true) {
    const ok = await new Promise<boolean>((resolve) => {
      const sock = createConnection({ host: HOST, port });
      const done = (v: boolean) => {
        sock.destroy();
        resolve(v);
      };
      sock.once("connect", () => done(true));
      sock.once("error", () => done(false));
    });
    if (ok) return;
    if (Date.now() > deadline) {
      throw new Error(
        `tauri-plugin-playwright control server not reachable on ${HOST}:${port} ` +
          `within ${timeoutMs}ms. Start the app with: bun run e2e:app`,
      );
    }
    await new Promise((r) => setTimeout(r, 500));
  }
}

export const test = base.extend<{ tauriPage: TauriPage }>({
  tauriPage: async ({}, use) => {
    await waitForPort(TCP_PORT);

    const client = new PluginClient(undefined, TCP_PORT);
    await client.connect();

    const ping = await client.send({ type: "ping" });
    if (!ping.ok) {
      client.disconnect();
      throw new Error(`plugin ping failed: ${ping.error ?? "no response"}`);
    }

    const page = new TauriPage(client);

    // Wait until the webview has finished loading and the plugin's init script
    // has run (it sets window.__PW_ACTIVE__ on every navigation).
    await page.waitForFunction(
      'document.readyState === "complete" && !!window.__PW_ACTIVE__',
      30_000,
    );

    try {
      await use(page);
    } finally {
      client.disconnect();
    }
  },
});

export const expect = tauriExpect;
