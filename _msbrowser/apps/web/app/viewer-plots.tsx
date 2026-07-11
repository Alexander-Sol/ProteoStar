"use client";

import dynamic from "next/dynamic";

import type { SpectrumPlotProps, TicPlotProps } from "@msbrowser/plot-adapter";

export const ClientTicPlot = dynamic<TicPlotProps>(
  async () => {
    const module = await import("@msbrowser/plot-adapter");
    return module.TicPlot;
  },
  {
    ssr: false
  }
);

export const ClientSpectrumPlot = dynamic<SpectrumPlotProps>(
  async () => {
    const module = await import("@msbrowser/plot-adapter");
    return module.SpectrumPlot;
  },
  {
    ssr: false
  }
);
