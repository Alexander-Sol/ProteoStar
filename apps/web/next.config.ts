import type { NextConfig } from "next";

const nextConfig: NextConfig = {
  transpilePackages: [
    "@msbrowser/imsp-core",
    "@msbrowser/viewer-state",
    "@msbrowser/plot-adapter",
    "@msbrowser/ui"
  ]
};

export default nextConfig;
