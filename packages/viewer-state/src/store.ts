import { createStore, type StoreApi } from "zustand/vanilla";

import { reduceViewerState } from "./reducer";
import { createInitialViewerState } from "./state";
import type { ViewerCommand, ViewerStoreState } from "./types";

export type ViewerStore = StoreApi<ViewerStoreState>;

export function createViewerStore(): ViewerStore {
  return createStore<ViewerStoreState>()((set) => ({
    ...createInitialViewerState(),
    dispatch(command: ViewerCommand) {
      set((state) => reduceViewerState(state, command));
    }
  }));
}
