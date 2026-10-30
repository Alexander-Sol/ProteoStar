import { describe, expect, it } from "vitest";

import {
  MODIFICATIONS,
  cleanSequence,
  fullSequence,
  predictMeanCharge,
  invalidResidues,
  modsFor,
  monoisotopicMass
} from "./proteoform-mass";

const UBIQUITIN =
  "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG";
const mod = (name: string) => MODIFICATIONS.find((m) => m.name === name)!;

describe("monoisotopicMass", () => {
  it("matches ubiquitin's published monoisotopic mass", () => {
    expect(monoisotopicMass(UBIQUITIN, [])).toBeCloseTo(8559.6167, 3);
  });

  it("adds modifications and disulfides", () => {
    const base = monoisotopicMass("ACKC", []);
    const mods = [
      { position: "N-term" as const, mod: mod("Acetylation") },
      { position: 2, mod: mod("Trimethylation") }
    ];
    expect(monoisotopicMass("ACKC", mods, 1)).toBeCloseTo(base + 42.01056 + 42.04695 - 2.01565, 4);
  });

  it("is NaN for a sequence with a non-residue letter", () => {
    expect(monoisotopicMass("PEPTIDEB", [])).toBeNaN();
  });
});

describe("predictMeanCharge", () => {
  it("is mzLib's ResidueChargeModel", () => {
    // Ubiquitin: 7 K, 4 R, 1 H, 5 D, 6 E.
    const m = monoisotopicMass(UBIQUITIN, []);
    const expected =
      -13.146636130364445 +
      0.22788186941351907 * Math.sqrt(m) +
      7 * 0.1463542063398008 +
      4 * 0.22587238463996417 +
      0.1808826550146898 -
      5 * 0.043768121160745424 -
      6 * 0.012663575728241415;
    expect(predictMeanCharge(UBIQUITIN, m)).toBeCloseTo(expected, 10);
    expect(predictMeanCharge(UBIQUITIN, m)).toBeCloseTo(9.75, 1);
  });
});

describe("sequence input", () => {
  it("strips a FASTA header, whitespace and numbering", () => {
    expect(cleanSequence(">sp|P0CG48|UBC\n  1 mqifv ktltg\n11 KTITL\r\n")).toBe("MQIFVKTLTGKTITL");
  });

  it("reports letters that are not residues", () => {
    expect(invalidResidues("PEPBXB")).toEqual(["B", "X"]);
  });
});

describe("modifications", () => {
  it("offers only the modifications a residue or terminus takes", () => {
    const names = (pos: number | "N-term") => modsFor("SKM", pos).map((m) => m.name);
    expect(names(0)).toEqual(["Phosphorylation"]);
    expect(names(2)).toEqual(["Oxidation"]);
    expect(names(1)).toContain("Acetylation");
    expect(names("N-term")).toEqual(["Acetylation", "Formylation", "Myristoylation"]);
  });

  it("writes the full sequence the way MetaMorpheus does", () => {
    const mods = [
      { position: "N-term" as const, mod: mod("Acetylation") },
      { position: 1, mod: mod("Methylation") }
    ];
    expect(fullSequence("SKM", mods)).toBe(
      "[Common Biological:Acetylation on X]SK[Common Biological:Methylation on K]M"
    );
  });
});
