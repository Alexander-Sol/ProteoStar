import { createImspDatasetProvider, type DatasetProvider } from "@msbrowser/imsp-core";

import type { ImspWorkerRequest, ImspWorkerResponse } from "./imsp-worker-protocol";

let provider: DatasetProvider | null = null;

self.onmessage = async (event: MessageEvent<ImspWorkerRequest>) => {
  const message = event.data;

  try {
    switch (message.type) {
      case "initialize": {
        provider = createImspDatasetProvider(message.buffer);
        postWorkerMessage({ id: message.id, ok: true, type: "initialize", result: true });
        return;
      }

      case "getMetadata":
        postWorkerMessage({
          id: message.id,
          ok: true,
          type: "getMetadata",
          result: await requireProvider().getMetadata()
        });
        return;

      case "getScanSummaries":
        postWorkerMessage({
          id: message.id,
          ok: true,
          type: "getScanSummaries",
          result: await requireProvider().getScanSummaries()
        });
        return;

      case "getNearestScan":
        postWorkerMessage({
          id: message.id,
          ok: true,
          type: "getNearestScan",
          result: await requireProvider().getNearestScan(message.retentionTime)
        });
        return;

      case "getSpectrumForScan":
        postWorkerMessage({
          id: message.id,
          ok: true,
          type: "getSpectrumForScan",
          result: await requireProvider().getSpectrumForScan(message.scanIndex)
        });
        return;

      case "getPeaksInMzRange":
        postWorkerMessage({
          id: message.id,
          ok: true,
          type: "getPeaksInMzRange",
          result: await requireProvider().getPeaksInMzRange(message.mzMin, message.mzMax)
        });
        return;

      case "getTicTrace":
        postWorkerMessage({
          id: message.id,
          ok: true,
          type: "getTicTrace",
          result: await requireProvider().getTicTrace(message.rtRange)
        });
        return;
    }
  } catch (error) {
    postWorkerMessage({
      id: message.id,
      ok: false,
      error: {
        name: error instanceof Error ? error.name : "Error",
        message: error instanceof Error ? error.message : "Worker request failed"
      }
    });
  }
};

function requireProvider(): DatasetProvider {
  if (!provider) {
    throw new Error("IMSP worker provider has not been initialized");
  }

  return provider;
}

function postWorkerMessage(message: ImspWorkerResponse): void {
  self.postMessage(message);
}
