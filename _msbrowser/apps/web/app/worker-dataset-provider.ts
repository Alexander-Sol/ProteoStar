import { createImspDatasetProvider, type DatasetProvider } from "@msbrowser/imsp-core";

import type { ImspWorkerRequest, ImspWorkerResponse } from "./imsp-worker-protocol";

type WorkerRequestWithoutId = ImspWorkerRequest extends infer Request
  ? Request extends { id: number }
    ? Omit<Request, "id">
    : never
  : never;

export interface BrowserDatasetProviderHandle {
  provider: DatasetProvider;
  dispose?: () => void;
}

export async function createBrowserDatasetProvider(
  buffer: ArrayBuffer
): Promise<BrowserDatasetProviderHandle> {
  if (typeof Worker === "undefined") {
    return {
      provider: createImspDatasetProvider(buffer)
    };
  }

  try {
    const worker = new Worker(new URL("./imsp-worker.ts", import.meta.url), {
      type: "module"
    });
    const rpc = createWorkerRpc(worker);

    await rpc.send({
      type: "initialize",
      buffer
    });

    return {
      provider: {
        getMetadata: () => rpc.send({ type: "getMetadata" }),
        getScanSummaries: () => rpc.send({ type: "getScanSummaries" }),
        getNearestScan: (retentionTime) =>
          rpc.send({ type: "getNearestScan", retentionTime }),
        getSpectrumForScan: (scanIndex) =>
          rpc.send({ type: "getSpectrumForScan", scanIndex }),
        getPeaksInMzRange: (mzMin, mzMax) =>
          rpc.send({ type: "getPeaksInMzRange", mzMin, mzMax }),
        getTicTrace: (rtRange) => rpc.send({ type: "getTicTrace", rtRange })
      },
      dispose: () => {
        rpc.dispose();
      }
    };
  } catch {
    return {
      provider: createImspDatasetProvider(buffer)
    };
  }
}

function createWorkerRpc(worker: Worker) {
  let nextId = 1;
  const pending = new Map<
    number,
    {
      resolve(value: unknown): void;
      reject(reason: unknown): void;
    }
  >();

  worker.onmessage = (event: MessageEvent<ImspWorkerResponse>) => {
    const message = event.data;
    const deferred = pending.get(message.id);

    if (!deferred) {
      return;
    }

    pending.delete(message.id);

    if (message.ok) {
      deferred.resolve(message.result);
      return;
    }

    deferred.reject(new Error(message.error.message));
  };

  worker.onerror = (event) => {
    const error = event.error instanceof Error ? event.error : new Error(event.message);
    for (const deferred of pending.values()) {
      deferred.reject(error);
    }
    pending.clear();
  };

  return {
    send<T>(
      request: WorkerRequestWithoutId,
      transfer: Transferable[] = []
    ): Promise<T> {
      return new Promise<T>((resolve, reject) => {
        const id = nextId;
        nextId += 1;
        pending.set(id, { resolve, reject });
        worker.postMessage({ ...request, id } as ImspWorkerRequest, transfer);
      });
    },
    dispose(): void {
      worker.terminate();
      for (const deferred of pending.values()) {
        deferred.reject(new Error("IMSP worker provider was terminated"));
      }
      pending.clear();
    }
  };
}
