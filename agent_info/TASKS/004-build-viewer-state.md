# Task 004: Build Viewer State

## Goal

Implement the state model for dataset loading, selection, pinning, and zoom behavior.

## Scope

1. Add dataset loading state.
2. Add TIC viewport and pin state.
3. Add spectrum viewport and pin state.
4. Add semantic commands for selection and reset.
5. Add tests for the core interaction rules.

## Deliverables

1. Zustand store or equivalent state module
2. State-driven interaction rules independent of Plotly
3. Unit tests for pin and zoom logic

## Notes

The store should consume semantic plot events, not raw Plotly event payloads.
