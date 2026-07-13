// Feature-finding walkthrough: request the top-down charge-state-ladder fit for a
// manually chosen seed and anchoring charge (`score_seed_ladder` command).

import { invoke } from "@tauri-apps/api/core";

import type { SeedLadder } from "./contract";

/**
 * Score one seed's full charge-state ladder for a given anchoring charge `zSeed`.
 * The backend takes a retention time (not a scan index) and resolves the nearest
 * MS1 scan in its own index space, so the seed lookup can't drift from the
 * displayed spectrum. Requires the peak index to be built.
 */
export function scoreSeedLadder(
  handle: number,
  retentionTime: number,
  seedMz: number,
  zSeed: number,
  maxCharge = 30
): Promise<SeedLadder> {
  return invoke<SeedLadder>("score_seed_ladder", {
    handle,
    retentionTime,
    seedMz,
    zSeed,
    options: { maxCharge }
  });
}
