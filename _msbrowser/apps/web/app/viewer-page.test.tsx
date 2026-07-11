import { readFileSync } from "node:fs";
import { join } from "node:path";

import { cleanup, fireEvent, render, screen, waitFor } from "@testing-library/react";
import React from "react";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";

vi.mock("@msbrowser/plot-adapter", () => ({}));

vi.mock("./viewer-plots", () => ({
  ClientTicPlot: ({ onEvent, traces, rangeSelectionEnabled, viewport }: any) => (
    <div>
      <div data-testid="tic-selected-0">{traces[0]?.selectedScanIndex ?? "none"}</div>
      <div data-testid="tic-range-enabled">{String(rangeSelectionEnabled)}</div>
      <div data-testid="tic-viewport">
        {viewport.xMin ?? "null"}:{viewport.xMax ?? "null"}
      </div>
      <button
        onClick={() =>
          onEvent({ type: "area-click", retentionTime: traces[0]?.points[2]?.retentionTime ?? 3 })
        }
        type="button"
      >
        Select Third TIC Point
      </button>
      <button
        onClick={() => onEvent({ type: "range-select", range: { min: 1.1, max: 2.2 } })}
        type="button"
      >
        Zoom TIC Range
      </button>
    </div>
  ),
  ClientSpectrumPlot: ({ onEvent, traces, rangeSelectionEnabled, viewport }: any) => (
    <div>
      <div data-testid="spectrum-peak-count">{traces[0]?.peaks?.length ?? 0}</div>
      <div data-testid="spectrum-range-enabled">{String(rangeSelectionEnabled)}</div>
      <div data-testid="spectrum-viewport">
        {viewport.xMin ?? "null"}:{viewport.xMax ?? "null"}
      </div>
      <button
        onClick={() => onEvent({ type: "range-select", range: { min: 700, max: 900 } })}
        type="button"
      >
        Zoom Spectrum Range
      </button>
    </div>
  )
}));

const { ViewerPage } = await import("./viewer-page");

describe("ViewerPage", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
  });

  afterEach(() => {
    cleanup();
  });

  it("loads a local IMSP file and renders the initial selected spectrum", async () => {
    render(<ViewerPage />);

    await uploadFixture("tiny-known.imsp");

    expect(await screen.findByText("tiny-known.imsp")).toBeTruthy();
    expect(await screen.findByText("1.000 min")).toBeTruthy();
    expect(screen.getByText("1")).toBeTruthy();
    expect((await screen.findByTestId("spectrum-peak-count")).textContent).toBe("3");
    expect(screen.getByTestId("tic-selected-0").textContent).toBe("0");
    // TIC is pinned by default; spectrum is not.
    expect(screen.getByTestId("tic-range-enabled").textContent).toBe("true");
    expect(screen.getByTestId("spectrum-range-enabled").textContent).toBe("false");
    expect(screen.getByTestId("spectrum-viewport").textContent).toBe("200:1200");
    expect(screen.getByText("200.0000 to 1200.0000")).toBeTruthy();
  });

  it("applies TIC zoom immediately since TIC is pinned by default", async () => {
    render(<ViewerPage />);

    await uploadFixture("tiny-known.imsp");
    await screen.findByTestId("tic-selected-0");

    fireEvent.click(screen.getByRole("button", { name: "Zoom TIC Range" }));

    await waitFor(() => {
      expect(screen.getByText("1.100 min to 2.200 min")).toBeTruthy();
    });
  });

  it("applies spectrum zoom only after the spectrum panel is pinned", async () => {
    render(<ViewerPage />);

    await uploadFixture("tiny-known.imsp");
    await screen.findByTestId("tic-selected-0");

    // Spectrum is not pinned — zoom should be ignored.
    fireEvent.click(screen.getByRole("button", { name: "Zoom Spectrum Range" }));
    expect(screen.getAllByText("Full range").length).toBeGreaterThan(0);

    // Pin the spectrum panel (second "Pin Zoom" button).
    const buttons = screen.getAllByRole("button", { name: "Pin Zoom" });
    fireEvent.click(buttons[1]);
    fireEvent.click(screen.getByRole("button", { name: "Zoom Spectrum Range" }));

    await waitFor(() => {
      expect(screen.getByText("700.0000 to 900.0000")).toBeTruthy();
    });
  });

  it("selecting a TIC scan updates the selected spectrum", async () => {
    render(<ViewerPage />);

    await uploadFixture("tiny-known.imsp");
    await screen.findByTestId("tic-selected-0");
    await screen.findByTestId("spectrum-peak-count");

    fireEvent.click(screen.getByRole("button", { name: "Select Third TIC Point" }));

    await waitFor(() => {
      expect(screen.getByTestId("tic-selected-0").textContent).toBe("2");
      expect(screen.getByText("3.000 min")).toBeTruthy();
    });
  });
});

async function uploadFixture(name: string): Promise<void> {
  const input = screen.getByLabelText("Open `.imsp` File") as HTMLInputElement;
  const filePath = join(process.cwd(), "fixtures", name);
  const bytes = readFileSync(filePath);
  const file = new File([bytes], name, { type: "application/octet-stream" });
  Object.defineProperty(file, "arrayBuffer", {
    value: async () => bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength)
  });

  fireEvent.change(input, {
    target: {
      files: [file]
    }
  });

  await waitFor(() => {
    expect(input.files?.[0]?.name).toBe(name);
  });
}
