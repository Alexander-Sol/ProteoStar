import { defineConfig } from "vitest/config";

export default defineConfig({
  test: {
    include: ["packages/*/tests/**/*.test.ts", "apps/desktop/**/*.test.tsx"],
    passWithNoTests: true,
    environment: "node",
    environmentMatchGlobs: [["apps/desktop/**/*.test.tsx", "jsdom"]]
  }
});
