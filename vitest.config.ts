import { defineConfig } from "vitest/config";

export default defineConfig({
  test: {
    include: ["packages/*/tests/**/*.test.ts", "apps/web/**/*.test.tsx"],
    environment: "node",
    environmentMatchGlobs: [["apps/web/**/*.test.tsx", "jsdom"]]
  }
});
