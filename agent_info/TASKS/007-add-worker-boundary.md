# Task 007: Add Worker Boundary

## Goal

Move parsing and expensive queries off the main UI thread.

## Scope

1. Create a worker entry point.
2. Move parse and selected query operations into the worker.
3. Define request and response messages.
4. Keep the UI API stable while switching execution context.

## Deliverables

1. Worker-backed dataset provider
2. Loading and error handling for async operations
3. No major UI freeze during medium-file load

## Notes

Keep the dataset provider interface stable so the UI does not care whether data is local-sync, local-worker, or backend-served.
