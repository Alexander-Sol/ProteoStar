import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// Vite config tuned for a Tauri desktop front end.
// - fixed dev port so Tauri's devUrl matches,
// - clearScreen off so cargo/tauri logs stay visible.
export default defineConfig({
  plugins: [react()],
  clearScreen: false,
  server: {
    port: 1420,
    strictPort: true
  }
});
