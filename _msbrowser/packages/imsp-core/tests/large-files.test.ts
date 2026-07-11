import { existsSync, readdirSync, readFileSync } from "node:fs";
import { join } from "node:path";
import { describe, expect, it } from "vitest";

import { createImspDatasetProvider } from "../src/index";

const largeFilesDirectory = join(process.cwd(), "LargeFiles");

function loadLargeFile(name: string): ArrayBuffer {
  const filePath = join(largeFilesDirectory, name);
  const file = readFileSync(filePath);
  return file.buffer.slice(file.byteOffset, file.byteOffset + file.byteLength);
}

describe("large IMSP file smoke coverage", () => {
  const largeFiles = existsSync(largeFilesDirectory)
    ? readdirSync(largeFilesDirectory).filter((fileName) => fileName.endsWith(".imsp"))
    : [];

  for (const fileName of largeFiles) {
    it(`loads ${fileName} metadata and TIC trace`, async () => {
      const dataset = createImspDatasetProvider(loadLargeFile(fileName));
      const [metadata, ticTrace] = await Promise.all([
        dataset.getMetadata(),
        dataset.getTicTrace()
      ]);

      expect(metadata.scanCount).toBeGreaterThan(0);
      expect(metadata.totalPeakCount).toBeGreaterThan(0);
      expect(metadata.mzRange).not.toBeNull();
      expect(ticTrace.length).toBe(metadata.scanCount);
    });
  }
});
