// Monoisotopic mass of a proteoform written as a sequence plus modifications: what the workbench's
// "Add proteoform" dialog computes to place a new proteoform in the simulation.

/** Monoisotopic residue masses (Da): the amino acid less one water. */
export const RESIDUE_MASS: Readonly<Record<string, number>> = {
  G: 57.02146372,
  A: 71.03711381,
  S: 87.0320284,
  P: 97.05276388,
  V: 99.06841395,
  T: 101.04767846,
  C: 103.00918451,
  L: 113.08406402,
  I: 113.08406402,
  N: 114.04292744,
  D: 115.02694303,
  Q: 128.0585775,
  K: 128.09496302,
  E: 129.04259309,
  M: 131.04048508,
  H: 137.05891186,
  F: 147.06841391,
  U: 150.95363559,
  R: 156.10111103,
  Y: 163.06332854,
  W: 186.07931295,
  O: 237.14772496
};

export const WATER_MASS = 18.0105646863;
/** One disulfide bond: two cysteine thiols lose a hydrogen each. */
export const DISULFIDE_MASS = -2.01565006;

/** Where a modification may sit. */
export type ModSite = { residues: string } | "N-term" | "C-term";

export interface ModDefinition {
  /** MetaMorpheus-style category, e.g. "Common Biological". */
  category: string;
  name: string;
  /** Monoisotopic mass change, Da. */
  delta: number;
  sites: readonly ModSite[];
}

export const MODIFICATIONS: readonly ModDefinition[] = [
  { category: "Common Biological", name: "Acetylation", delta: 42.01056469, sites: [{ residues: "K" }, "N-term"] },
  { category: "Common Biological", name: "Methylation", delta: 14.01565006, sites: [{ residues: "KR" }] },
  { category: "Common Biological", name: "Dimethylation", delta: 28.03130013, sites: [{ residues: "KR" }] },
  { category: "Common Biological", name: "Trimethylation", delta: 42.04695019, sites: [{ residues: "K" }] },
  { category: "Common Biological", name: "Phosphorylation", delta: 79.96633052, sites: [{ residues: "STY" }] },
  { category: "Common Biological", name: "Succinylation", delta: 100.01604399, sites: [{ residues: "K" }] },
  { category: "Common Biological", name: "Crotonylation", delta: 68.02621475, sites: [{ residues: "K" }] },
  { category: "Common Biological", name: "Formylation", delta: 27.99491462, sites: [{ residues: "K" }, "N-term"] },
  { category: "Common Biological", name: "Hydroxylation", delta: 15.99491462, sites: [{ residues: "KP" }] },
  { category: "Common Biological", name: "Citrullination", delta: 0.98401558, sites: [{ residues: "R" }] },
  { category: "Common Biological", name: "GlyGly", delta: 114.04292744, sites: [{ residues: "K" }] },
  { category: "Common Biological", name: "Myristoylation", delta: 210.19837021, sites: [{ residues: "G" }, "N-term"] },
  { category: "Common Biological", name: "Palmitoylation", delta: 238.22966834, sites: [{ residues: "C" }] },
  { category: "Common Biological", name: "Amidation", delta: -0.98401558, sites: ["C-term"] },
  { category: "Common Artifact", name: "Deamidation", delta: 0.98401558, sites: [{ residues: "NQ" }] },
  { category: "Common Artifact", name: "Pyro-glu from Q", delta: -17.02654911, sites: [{ residues: "Q" }] },
  { category: "Common Artifact", name: "Pyro-glu from E", delta: -18.01056469, sites: [{ residues: "E" }] },
  { category: "Common Variable", name: "Oxidation", delta: 15.99491462, sites: [{ residues: "M" }] },
  { category: "Common Fixed", name: "Carbamidomethyl", delta: 57.02146372, sites: [{ residues: "C" }] }
];

/** Where a placed modification sits: a 0-based residue index, or a terminus. */
export type ModPosition = number | "N-term" | "C-term";

export interface PlacedMod {
  position: ModPosition;
  mod: ModDefinition;
}

/** The modifications that may sit at `position` of `sequence`. */
export function modsFor(sequence: string, position: ModPosition): ModDefinition[] {
  return MODIFICATIONS.filter((m) =>
    m.sites.some((s) =>
      typeof s === "string"
        ? s === position
        : typeof position === "number" && s.residues.includes(sequence[position] ?? "")
    )
  );
}

/**
 * A pasted sequence reduced to its residue letters: whitespace, digits and FASTA line breaks are
 * dropped and letters upper-cased. A leading FASTA header line is dropped too.
 */
export function cleanSequence(text: string): string {
  const lines = text.split(/\r?\n/).filter((l) => !l.trimStart().startsWith(">"));
  return lines.join("").replace(/[\s\d]/g, "").toUpperCase();
}

/** Characters in a cleaned sequence that are not residues, each once. */
export function invalidResidues(sequence: string): string[] {
  return [...new Set([...sequence].filter((c) => RESIDUE_MASS[c] === undefined))];
}

/** Monoisotopic neutral mass: residues + water + modifications + disulfide bonds. */
export function monoisotopicMass(sequence: string, mods: readonly PlacedMod[], disulfides = 0): number {
  let mass = WATER_MASS;
  for (const c of sequence) mass += RESIDUE_MASS[c] ?? NaN;
  for (const m of mods) mass += m.mod.delta;
  return mass + disulfides * DISULFIDE_MASS;
}

/**
 * mzLib's `ResidueChargeModel` (TopDownSimulator/Prediction/IdOnlyPriors.cs): the centre of a
 * proteoform's charge envelope from √mass and a weight per charged residue, counted in the
 * unmodified sequence. Ordinary least squares on 309 Jurkat rep2 fract7 species fitted over their
 * whole charge range (`idonly-priors…rep2_fract7.v5.json`, residual SD 1.52 charges); on 329
 * held-out rep1 fract7 proteoforms, R² 0.70 and mean |error| 0.88 charges.
 */
export const CHARGE_MODEL = {
  intercept: -13.146636130364445,
  sqrtMassSlope: 0.22788186941351907,
  residueWeights: {
    K: 0.1463542063398008,
    R: 0.22587238463996417,
    H: 0.1808826550146898,
    D: -0.043768121160745424,
    E: -0.012663575728241415
  } as Readonly<Record<string, number>>
} as const;

export const CHARGE_MODEL_DESCRIPTION =
  "Linear model (mzLib ResidueChargeModel): mean charge = −13.15 + 0.228·√(mass in Da) " +
  "+ 0.146·K + 0.226·R + 0.181·H − 0.044·D − 0.013·E, counting each residue in the unmodified " +
  "sequence. Arginine and histidine add more charge than lysine; the acidic residues barely " +
  "matter. Trained on 309 Jurkat rep2 fraction 7 proteoforms; on 329 held-out rep1 fraction 7 " +
  "proteoforms R² 0.70, mean |error| 0.88 charges.";

/** {@link CHARGE_MODEL}'s mean charge for a proteoform of this sequence and mass. */
export function predictMeanCharge(sequence: string, massDa: number): number {
  const { intercept, sqrtMassSlope, residueWeights } = CHARGE_MODEL;
  let mu = intercept + sqrtMassSlope * Math.sqrt(massDa);
  for (const c of sequence) mu += residueWeights[c] ?? 0;
  return mu;
}

/**
 * The proteoform written the way MetaMorpheus writes a full sequence: each modification in
 * brackets after its residue (an N-terminal one before the first residue), e.g.
 * `[Common Biological:Acetylation on X]SGRGK[Common Biological:Methylation on K]…`.
 */
export function fullSequence(sequence: string, mods: readonly PlacedMod[]): string {
  const tag = (m: PlacedMod, motif: string) => `[${m.mod.category}:${m.mod.name} on ${motif}]`;
  const nTerm = mods.filter((m) => m.position === "N-term").map((m) => tag(m, "X"));
  const cTerm = mods.filter((m) => m.position === "C-term").map((m) => tag(m, "X"));
  let out = nTerm.join("");
  for (let i = 0; i < sequence.length; i++) {
    out += sequence[i];
    for (const m of mods) if (m.position === i) out += tag(m, sequence[i]);
  }
  return out + (cTerm.length > 0 ? `-${cTerm.join("")}` : "");
}
