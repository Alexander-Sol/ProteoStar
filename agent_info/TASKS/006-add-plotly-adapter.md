# Task 006: Add Plotly Adapter

## Goal

Wrap Plotly behind a small adapter layer so plotting can be replaced later if needed.

## Scope

1. Define plot adapter interfaces.
2. Implement TIC rendering through Plotly.
3. Implement spectrum rendering through Plotly.
4. Translate Plotly interactions into semantic app events.

## Deliverables

1. Plot adapter interface
2. Plotly-backed adapter implementation
3. No Plotly event leakage into viewer-state logic

## Notes

This is the main hedge against future migration to `uPlot` or a more specialized renderer.
