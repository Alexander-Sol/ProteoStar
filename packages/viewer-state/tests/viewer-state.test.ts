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

describe("viewer state reducer", () => {
  it("initializes dataset loading state", () => {
    expect(createInitialViewerState()).toEqual({
      activePanel: "tic",
      dataset: {
        status: "idle",
        dataset: null,
        errorMessage: null
      },
      selectedScanIndex: null,
      ticPanel: { pinned: false, range: null },
      spectrumPanel: { pinned: true, range: null }
    });
  });

  it("resets selection and panel state when dataset load succeeds", () => {
    const loadingState = reduceViewerState(createInitialViewerState(), {
      type: "dataset/load-started"
    });
    const readyState = reduceViewerState(loadingState, {
      type: "dataset/load-succeeded",
      dataset
    });

    expect(readyState.dataset).toEqual({
      status: "ready",
      dataset,
      errorMessage: null
    });
    expect(readyState.selectedScanIndex).toBeNull();
    expect(readyState.ticPanel).toEqual({ pinned: false, range: null });
    expect(readyState.spectrumPanel).toEqual({ pinned: true, range: null });
  });

  it("selects the nearest scan without mutating panel view state", () => {
    const state = reduceViewerState(
      reduceViewerState(createInitialViewerState(), {
        type: "dataset/load-succeeded",
        dataset
      }),
      { type: "panel/set-pinned", panelId: "tic", pinned: true }
    );

    const nextState = reduceViewerState(state, {
      type: "selection/select-nearest-scan",
      retentionTime: 1.51
    });

    expect(nextState.selectedScanIndex).toBe(1);
    expect(nextState.ticPanel).toEqual({ pinned: true, range: null });
    expect(nextState.spectrumPanel).toEqual({ pinned: true, range: null });
  });

  it("only zooms a pinned TIC panel", () => {
    const readyState = reduceViewerState(createInitialViewerState(), {
      type: "dataset/load-succeeded",
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
      dataset
    });
    state = reduceViewerState(state, { type: "selection/set-scan", scanIndex: 2 });
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

    expect(nextState.selectedScanIndex).toBe(2);
    expect(nextState.ticPanel).toEqual({ pinned: true, range: null });
    expect(nextState.spectrumPanel).toEqual({
      pinned: true,
      range: { min: 700, max: 900 }
    });
  });

  it("throws for scan selection before dataset load or out-of-range selection", () => {
    expect(() =>
      reduceViewerState(createInitialViewerState(), {
        type: "selection/select-nearest-scan",
        retentionTime: 1.5
      })
    ).toThrowError(ViewerStateError);

    const readyState = reduceViewerState(createInitialViewerState(), {
      type: "dataset/load-succeeded",
      dataset
    });

    expect(() =>
      reduceViewerState(readyState, { type: "selection/set-scan", scanIndex: 3 })
    ).toThrowError(ViewerStateError);
  });
});

describe("viewer store", () => {
  it("dispatches commands through the zustand-backed store", () => {
    const store = createViewerStore();

    store.getState().dispatch({ type: "dataset/load-succeeded", dataset });
    store.getState().dispatch({
      type: "selection/select-nearest-scan",
      retentionTime: 2.6
    });

    expect(store.getState().selectedScanIndex).toBe(2);
  });
});
