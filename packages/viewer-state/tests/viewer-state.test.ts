import { describe, expect, it } from "vitest";

import {
  ViewerStateError,
  createInitialViewerState,
  createViewerStore,
  reduceViewerState,
  type ViewerDataset
} from "../src/index";

const dataset: ViewerDataset = {
  metadata: {
    retentionTimeRange: { min: 1, max: 3 },
    mzRange: { min: 500.1234, max: 1250.3456 },
    scanCount: 3
  },
  scanSummaries: [
    { scanIndex: 0, oneBasedScanNumber: 1, retentionTime: 1, tic: 225000 },
    { scanIndex: 1, oneBasedScanNumber: 2, retentionTime: 2, tic: 250000 },
    { scanIndex: 2, oneBasedScanNumber: 3, retentionTime: 3, tic: 210000 }
  ]
};

const idleSlot = {
  load: { status: "idle", dataset: null, errorMessage: null },
  selectedScanIndex: null
};

describe("viewer state reducer", () => {
  it("initializes dataset loading state", () => {
    expect(createInitialViewerState()).toEqual({
      activePanel: "tic",
      datasetSlots: [idleSlot, idleSlot],
      ticPanel: { pinned: false, range: null },
      spectrumPanel: { pinned: true, range: null }
    });
  });

  it("resets selection and slot load state when dataset load succeeds", () => {
    const loadingState = reduceViewerState(createInitialViewerState(), {
      type: "dataset/load-started",
      slotIndex: 0
    });
    const readyState = reduceViewerState(loadingState, {
      type: "dataset/load-succeeded",
      slotIndex: 0,
      dataset
    });

    expect(readyState.datasetSlots[0].load).toEqual({
      status: "ready",
      dataset,
      errorMessage: null
    });
    expect(readyState.datasetSlots[0].selectedScanIndex).toBeNull();
    expect(readyState.datasetSlots[1]).toEqual(idleSlot);
  });

  it("resets panels only when slot 0 loads and slot 1 is idle", () => {
    // Load slot 0 first — should reset panels.
    const afterLoad0 = reduceViewerState(
      {
        ...createInitialViewerState(),
        ticPanel: { pinned: true, range: { min: 1, max: 2 } },
        spectrumPanel: { pinned: false, range: { min: 700, max: 900 } }
      },
      { type: "dataset/load-started", slotIndex: 0 }
    );
    expect(afterLoad0.ticPanel).toEqual({ pinned: false, range: null });
    expect(afterLoad0.spectrumPanel).toEqual({ pinned: true, range: null });

    // Load slot 1 while slot 0 is already loaded — panels must not change.
    const withSlot0Ready = reduceViewerState(afterLoad0, {
      type: "dataset/load-succeeded",
      slotIndex: 0,
      dataset
    });
    const zoomed = {
      ...withSlot0Ready,
      ticPanel: { pinned: true, range: { min: 1, max: 2 } }
    };
    const afterLoad1 = reduceViewerState(zoomed, {
      type: "dataset/load-started",
      slotIndex: 1
    });
    expect(afterLoad1.ticPanel).toEqual({ pinned: true, range: { min: 1, max: 2 } });
  });

  it("selects the nearest scan without mutating panel view state", () => {
    const state = reduceViewerState(
      reduceViewerState(createInitialViewerState(), {
        type: "dataset/load-succeeded",
        slotIndex: 0,
        dataset
      }),
      { type: "panel/set-pinned", panelId: "tic", pinned: true }
    );

    const nextState = reduceViewerState(state, {
      type: "selection/select-nearest-scan",
      slotIndex: 0,
      retentionTime: 1.51
    });

    expect(nextState.datasetSlots[0].selectedScanIndex).toBe(1);
    expect(nextState.ticPanel).toEqual({ pinned: true, range: null });
    expect(nextState.spectrumPanel).toEqual({ pinned: false, range: null });
  });

  it("only zooms a pinned TIC panel", () => {
    const readyState = reduceViewerState(createInitialViewerState(), {
      type: "dataset/load-succeeded",
      slotIndex: 0,
      dataset
    });

    const ignoredZoomState = reduceViewerState(readyState, {
      type: "panel/zoom",
      panelId: "tic",
      range: { min: 2.5, max: 1.5 }
    });

    expect(ignoredZoomState.ticPanel.range).toBeNull();

    const pinnedState = reduceViewerState(readyState, {
      type: "panel/set-pinned",
      panelId: "tic",
      pinned: true
    });

    const zoomedState = reduceViewerState(pinnedState, {
      type: "panel/zoom",
      panelId: "tic",
      range: { min: 2.5, max: 1.5 }
    });

    expect(zoomedState.ticPanel).toEqual({
      pinned: true,
      range: { min: 1.5, max: 2.5 }
    });
  });

  it("only zooms a pinned spectrum panel", () => {
    const readyState = reduceViewerState(createInitialViewerState(), {
      type: "dataset/load-succeeded",
      slotIndex: 0,
      dataset
    });

    const zoomedState = reduceViewerState(readyState, {
      type: "panel/zoom",
      panelId: "spectrum",
      range: { min: 1000.9, max: 750.56 }
    });

    expect(zoomedState.spectrumPanel).toEqual({
      pinned: true,
      range: { min: 750.56, max: 1000.9 }
    });
    expect(zoomedState.ticPanel.range).toBeNull();
  });

  it("resets a single panel without disturbing the other panel or selection", () => {
    let state = reduceViewerState(createInitialViewerState(), {
      type: "dataset/load-succeeded",
      slotIndex: 0,
      dataset
    });
    state = reduceViewerState(state, {
      type: "selection/set-scan",
      slotIndex: 0,
      scanIndex: 2
    });
    state = reduceViewerState(state, {
      type: "panel/set-pinned",
      panelId: "tic",
      pinned: true
    });
    state = reduceViewerState(state, {
      type: "panel/set-pinned",
      panelId: "spectrum",
      pinned: true
    });
    state = reduceViewerState(state, {
      type: "panel/zoom",
      panelId: "tic",
      range: { min: 1, max: 2 }
    });
    state = reduceViewerState(state, {
      type: "panel/zoom",
      panelId: "spectrum",
      range: { min: 700, max: 900 }
    });

    const nextState = reduceViewerState(state, {
      type: "panel/reset",
      panelId: "tic"
    });

    expect(nextState.datasetSlots[0].selectedScanIndex).toBe(2);
    expect(nextState.ticPanel).toEqual({ pinned: false, range: null });
    expect(nextState.spectrumPanel).toEqual({
      pinned: true,
      range: { min: 700, max: 900 }
    });
  });

  it("throws for scan selection before dataset load or out-of-range selection", () => {
    expect(() =>
      reduceViewerState(createInitialViewerState(), {
        type: "selection/select-nearest-scan",
        slotIndex: 0,
        retentionTime: 1.5
      })
    ).toThrowError(ViewerStateError);

    const readyState = reduceViewerState(createInitialViewerState(), {
      type: "dataset/load-succeeded",
      slotIndex: 0,
      dataset
    });

    expect(() =>
      reduceViewerState(readyState, {
        type: "selection/set-scan",
        slotIndex: 0,
        scanIndex: 3
      })
    ).toThrowError(ViewerStateError);
  });

  it("manages two dataset slots independently", () => {
    let state = createInitialViewerState();

    state = reduceViewerState(state, {
      type: "dataset/load-succeeded",
      slotIndex: 0,
      dataset
    });
    state = reduceViewerState(state, {
      type: "selection/set-scan",
      slotIndex: 0,
      scanIndex: 1
    });
    state = reduceViewerState(state, {
      type: "dataset/load-succeeded",
      slotIndex: 1,
      dataset
    });
    state = reduceViewerState(state, {
      type: "selection/set-scan",
      slotIndex: 1,
      scanIndex: 2
    });

    expect(state.datasetSlots[0].selectedScanIndex).toBe(1);
    expect(state.datasetSlots[1].selectedScanIndex).toBe(2);
    expect(state.datasetSlots[0].load.status).toBe("ready");
    expect(state.datasetSlots[1].load.status).toBe("ready");
  });
});

describe("viewer store", () => {
  it("dispatches commands through the zustand-backed store", () => {
    const store = createViewerStore();

    store.getState().dispatch({
      type: "dataset/load-succeeded",
      slotIndex: 0,
      dataset
    });
    store.getState().dispatch({
      type: "selection/select-nearest-scan",
      slotIndex: 0,
      retentionTime: 2.6
    });

    expect(store.getState().datasetSlots[0].selectedScanIndex).toBe(2);
  });
});
