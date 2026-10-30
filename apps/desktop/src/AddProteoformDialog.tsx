// "Add proteoform": a new simulated proteoform from its sequence and modifications. The mass comes
// from the sequence; the elution and charge shape start from the spec's own proteoforms (their
// medians), so a new one looks like its neighbours until it is tuned.

import { useEffect, useMemo, useState } from "react";
import type { CSSProperties } from "react";

import {
  CHARGE_MODEL_DESCRIPTION,
  cleanSequence,
  predictMeanCharge,
  fullSequence,
  invalidResidues,
  modsFor,
  monoisotopicMass,
  type ModPosition,
  type PlacedMod
} from "./proteoform-mass";
import type { SimProteoform, SimSpec } from "./sim-spec";
import { NumberField, SIM_COLOR, SIM_SELECTED_COLOR, SmallButton, TextField } from "./sim-ui";

/** Residues per row of the residue map, numbered at each row start. */
const RESIDUES_PER_ROW = 20;

function median(values: number[], fallback: number): number {
  const v = values.filter(Number.isFinite).sort((a, b) => a - b);
  return v.length === 0 ? fallback : v[v.length >> 1];
}

/** The shape a new proteoform starts with: the spec's medians, or generic top-down values. */
function templateOf(spec: SimSpec) {
  const ps = spec.proteoforms;
  return {
    abundance: median(ps.map((p) => p.abundance), 1e7),
    sigma: median(ps.map((p) => p.rt.sigma), 0.2),
    tau: median(ps.map((p) => p.rt.tau), 0.1),
    sigmaZ: median(ps.map((p) => p.charge.sigmaZ), 3),
    // Where the charge envelope centres: the typical mass-to-mean-charge ratio.
    mzAtMeanCharge: median(ps.map((p) => p.monoisotopicMass / p.charge.muZ), 800)
  };
}

const positionLabel = (sequence: string, pos: ModPosition): string =>
  typeof pos === "number" ? `${sequence[pos]}${pos + 1}` : pos;

export function AddProteoformDialog({
  spec,
  defaultRt,
  onAdd,
  onClose
}: {
  spec: SimSpec;
  defaultRt: number;
  onAdd: (p: SimProteoform) => void;
  onClose: () => void;
}) {
  const template = useMemo(() => templateOf(spec), [spec.proteoforms]);
  const [name, setName] = useState("");
  const [rawSequence, setRawSequence] = useState("");
  // Null follows the N-terminal methionine excision rule; set once the user ticks the box.
  const [removeMetChoice, setRemoveMet] = useState<boolean | null>(null);
  const [mods, setMods] = useState<PlacedMod[]>([]);
  const [disulfides, setDisulfides] = useState(0);
  const [picked, setPicked] = useState<ModPosition | null>(null);
  const [modChoice, setModChoice] = useState(0);
  const [abundance, setAbundance] = useState(template.abundance);
  const [rt, setRt] = useState(Number(defaultRt.toFixed(2)));
  const [rtSigma, setRtSigma] = useState(template.sigma);
  // Null follows the mass; set once the user moves it.
  const [meanCharge, setMeanCharge] = useState<number | null>(null);
  // What set the mean charge, when the sequence model did.
  const [chargeHint, setChargeHint] = useState<string | null>(null);

  const cleaned = cleanSequence(rawSequence);
  const startsWithMet = cleaned.startsWith("M") && cleaned.length > 1;
  // Methionine aminopeptidase cleaves the initiator Met when the next residue is small.
  const removeMet = removeMetChoice ?? "GASCTPV".includes(cleaned[1] ?? "");
  const sequence = startsWithMet && removeMet ? cleaned.slice(1) : cleaned;
  const invalid = invalidResidues(sequence);
  const cysteines = [...sequence].filter((c) => c === "C").length;
  const maxDisulfides = cysteines >> 1;

  // A sequence edit drops the modifications that no longer fit where they sit.
  useEffect(() => {
    setMods((ms) => ms.filter((m) => modsFor(sequence, m.position).some((d) => d === m.mod)));
    setPicked((p) => (typeof p === "number" && p >= sequence.length ? null : p));
    setDisulfides((d) => Math.min(d, maxDisulfides));
  }, [sequence, maxDisulfides]);

  const mass = sequence.length > 0 && invalid.length === 0 ? monoisotopicMass(sequence, mods, disulfides) : null;
  const { min: zMin, max: zMax } = spec.chargeRange;
  const autoCharge =
    mass !== null
      ? Math.min(zMax, Math.max(zMin, Math.round((2 * mass) / template.mzAtMeanCharge) / 2))
      : (zMin + zMax) / 2;
  const muZ = meanCharge ?? autoCharge;

  const choices = picked !== null ? modsFor(sequence, picked) : [];
  const modAt = (pos: ModPosition) => mods.find((m) => m.position === pos);
  const rtMin = spec.scanTimes[0] ?? 0;
  const rtMax = spec.scanTimes[spec.scanTimes.length - 1] ?? 60;

  const add = () => {
    if (mass === null) return;
    const full = fullSequence(sequence, mods);
    onAdd({
      monoisotopicMass: mass,
      abundance,
      rt: { mu: rt, sigma: rtSigma, tau: template.tau },
      charge: { kind: "Gaussian", muZ, sigmaZ: template.sigmaZ },
      identifier: name.trim() ? `${name.trim()}:${full}` : full
    });
  };

  // Escape closes, as a dialog should.
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if (e.key === "Escape") onClose();
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [onClose]);

  const terminus = (pos: "N-term" | "C-term") => (
    <button
      type="button"
      onClick={() => {
        setPicked(pos);
        setModChoice(0);
      }}
      title={modAt(pos) ? `${modAt(pos)!.mod.name}` : `${pos} — click to modify`}
      style={residueStyle(picked === pos, !!modAt(pos), true)}
    >
      {pos === "N-term" ? "N" : "C"}
    </button>
  );

  return (
    <div style={backdropStyle} onMouseDown={(e) => e.target === e.currentTarget && onClose()}>
      <div style={dialogStyle} role="dialog" aria-label="Add proteoform">
        <div style={titleRowStyle}>
          <strong>Add proteoform</strong>
          <button type="button" onClick={onClose} style={closeStyle} aria-label="Close">
            ✕
          </button>
        </div>

        <div style={bodyStyle}>
          <TextField label="Name (optional, e.g. an accession)" value={name} onChange={setName} />

          <label style={fieldStyle}>
            <span>Sequence</span>
            <textarea
              style={textareaStyle}
              rows={4}
              placeholder="Paste a sequence or a FASTA entry"
              value={rawSequence}
              onChange={(e) => setRawSequence(e.target.value)}
              spellCheck={false}
            />
          </label>
          <div style={noteStyle}>
            {cleaned.length === 0
              ? "Whitespace, numbering and a FASTA header are ignored."
              : invalid.length > 0
                ? <span style={errorStyle}>Not amino acids: {invalid.join(", ")}</span>
                : `${sequence.length} residues`}
          </div>
          {startsWithMet ? (
            <label style={checkStyle}>
              <input type="checkbox" checked={removeMet} onChange={(e) => setRemoveMet(e.target.checked)} />
              Remove the initiator methionine
              {removeMetChoice === null ? (
                <span style={noteStyle}>
                  ({removeMet ? "cleaved" : "kept"}: followed by {cleaned[1]})
                </span>
              ) : null}
            </label>
          ) : null}

          {sequence.length > 0 && invalid.length === 0 ? (
            <>
              <div style={sectionTitleStyle}>Modifications</div>
              <div style={noteStyle}>Click a residue or a terminus to modify it.</div>
              <div style={mapStyle}>
                {terminus("N-term")}
                {Array.from({ length: Math.ceil(sequence.length / RESIDUES_PER_ROW) }, (_, row) => (
                  <div key={row} style={mapRowStyle}>
                    <span style={rowNumberStyle}>{row * RESIDUES_PER_ROW + 1}</span>
                    {[...sequence.slice(row * RESIDUES_PER_ROW, (row + 1) * RESIDUES_PER_ROW)].map((c, k) => {
                      const i = row * RESIDUES_PER_ROW + k;
                      const m = modAt(i);
                      return (
                        <button
                          key={i}
                          type="button"
                          onClick={() => {
                            setPicked(i);
                            setModChoice(0);
                          }}
                          title={m ? `${c}${i + 1}: ${m.mod.name}` : `${c}${i + 1}`}
                          style={residueStyle(picked === i, !!m, false)}
                        >
                          {c}
                        </button>
                      );
                    })}
                  </div>
                ))}
                {terminus("C-term")}
              </div>

              {picked !== null ? (
                <div style={pickRowStyle}>
                  <span>{positionLabel(sequence, picked)}:</span>
                  {choices.length === 0 ? (
                    <span style={noteStyle}>no modifications defined for this residue</span>
                  ) : (
                    <>
                      <select
                        style={selectStyle}
                        value={modChoice}
                        onChange={(e) => setModChoice(Number(e.target.value))}
                      >
                        {choices.map((m, k) => (
                          <option key={m.name} value={k}>
                            {m.name} ({m.delta > 0 ? "+" : ""}
                            {m.delta.toFixed(4)})
                          </option>
                        ))}
                      </select>
                      <SmallButton
                        onClick={() => {
                          const mod = choices[modChoice];
                          if (!mod) return;
                          // One modification per site: a new one replaces the old.
                          setMods((ms) => [...ms.filter((m) => m.position !== picked), { position: picked, mod }]);
                        }}
                      >
                        {modAt(picked) ? "Replace" : "Add"}
                      </SmallButton>
                    </>
                  )}
                </div>
              ) : null}

              {mods.length > 0 ? (
                <ul style={modListStyle}>
                  {[...mods]
                    .sort((a, b) => order(a.position) - order(b.position))
                    .map((m) => (
                      <li key={String(m.position)} style={modItemStyle}>
                        <span>
                          {positionLabel(sequence, m.position)} · {m.mod.name}{" "}
                          <span style={noteStyle}>
                            ({m.mod.delta > 0 ? "+" : ""}
                            {m.mod.delta.toFixed(4)} Da)
                          </span>
                        </span>
                        <button
                          type="button"
                          style={closeStyle}
                          onClick={() => setMods((ms) => ms.filter((x) => x !== m))}
                          aria-label={`Remove ${m.mod.name}`}
                        >
                          ✕
                        </button>
                      </li>
                    ))}
                </ul>
              ) : null}

              {maxDisulfides > 0 ? (
                <NumberField
                  label="Disulfide bonds"
                  hint={`${cysteines} cysteines`}
                  value={disulfides}
                  step={1}
                  decimals={0}
                  onChange={(v) => setDisulfides(Math.max(0, Math.min(maxDisulfides, Math.round(v))))}
                />
              ) : null}
            </>
          ) : null}

          <div style={massStyle}>
            Monoisotopic mass:{" "}
            <strong>{mass !== null ? `${mass.toFixed(4)} Da` : "—"}</strong>
          </div>

          <div style={sectionTitleStyle}>In the simulation</div>
          <div style={pairStyle}>
            <NumberField
              label="Abundance"
              value={abundance}
              decimals={0}
              slider={{ min: 1e3, max: 1e12, log: true }}
              onChange={setAbundance}
            />
            <div>
              <NumberField
                label="Mean charge"
                hint={meanCharge === null ? "from the mass" : chargeHint ?? undefined}
                value={muZ}
                step={0.5}
                decimals={1}
                slider={{ min: zMin, max: zMax }}
                onChange={(v) => {
                  setMeanCharge(v);
                  setChargeHint(null);
                }}
              />
              <SmallButton
                onClick={() => {
                  if (mass === null) return;
                  const predicted = predictMeanCharge(sequence, mass);
                  setMeanCharge(Math.min(zMax, Math.max(zMin, Math.round(predicted * 2) / 2)));
                  // ± the model's residual SD on its training species.
                  setChargeHint(`predicted ${predicted.toFixed(1)} ± 1.5`);
                }}
                title={CHARGE_MODEL_DESCRIPTION}
              >
                Predict charge from sequence
              </SmallButton>
            </div>
            <NumberField
              label="Retention time (min)"
              hint="elution apex"
              value={rt}
              step={0.01}
              slider={{ min: rtMin, max: rtMax }}
              onChange={setRt}
            />
            <NumberField
              label="RT σ (min)"
              hint="Chromatographic Peak Width"
              value={rtSigma}
              step={0.01}
              slider={{ min: 0.01, max: 2 }}
              onChange={setRtSigma}
            />
          </div>
        </div>

        <div style={footerStyle}>
          <button type="button" style={cancelStyle} onClick={onClose}>
            Cancel
          </button>
          <button
            type="button"
            style={mass !== null ? addStyle : { ...addStyle, ...addDisabledStyle }}
            disabled={mass === null}
            onClick={add}
          >
            Add proteoform
          </button>
        </div>
      </div>
    </div>
  );
}

/** N-terminus first, then residues in order, then the C-terminus. */
const order = (pos: ModPosition): number => (pos === "N-term" ? -1 : pos === "C-term" ? Infinity : pos);

// ------------------------------------------------------------------ styles

const backdropStyle: CSSProperties = {
  position: "fixed",
  inset: 0,
  background: "rgba(20, 32, 48, 0.35)",
  display: "flex",
  alignItems: "center",
  justifyContent: "center",
  zIndex: 1000
};
const dialogStyle: CSSProperties = {
  width: "min(760px, 94vw)",
  maxHeight: "90vh",
  display: "flex",
  flexDirection: "column",
  background: "#ffffff",
  borderRadius: 8,
  boxShadow: "0 12px 40px rgba(20, 32, 48, 0.3)",
  color: "#24364d",
  fontSize: "0.82rem"
};
const titleRowStyle: CSSProperties = {
  display: "flex",
  justifyContent: "space-between",
  alignItems: "center",
  padding: "10px 16px",
  borderBottom: "1px solid #d5deeb",
  fontSize: "0.95rem"
};
const bodyStyle: CSSProperties = { overflow: "auto", padding: "8px 16px", minHeight: 0 };
const footerStyle: CSSProperties = {
  display: "flex",
  justifyContent: "flex-end",
  gap: 8,
  padding: "10px 16px",
  borderTop: "1px solid #d5deeb"
};
const closeStyle: CSSProperties = { border: "none", background: "transparent", cursor: "pointer", fontSize: "0.9rem" };
const fieldStyle: CSSProperties = { display: "flex", flexDirection: "column", gap: 2, margin: "6px 0", color: "#3a4a60" };
const textareaStyle: CSSProperties = {
  fontFamily: "ui-monospace, Consolas, monospace",
  fontSize: "0.8rem",
  padding: 6,
  border: "1px solid #cbd6e5",
  borderRadius: 4,
  resize: "vertical",
  wordBreak: "break-all"
};
const noteStyle: CSSProperties = { color: "#5a6b80", margin: "2px 0" };
const errorStyle: CSSProperties = { color: "#b02a37" };
const checkStyle: CSSProperties = { display: "flex", alignItems: "center", gap: 6, margin: "4px 0" };
const sectionTitleStyle: CSSProperties = {
  fontWeight: 600,
  marginTop: 12,
  paddingTop: 8,
  borderTop: "1px solid #e3e9f2"
};
const mapStyle: CSSProperties = {
  display: "flex",
  flexDirection: "column",
  alignItems: "flex-start",
  gap: 2,
  margin: "6px 0",
  fontFamily: "ui-monospace, Consolas, monospace"
};
const mapRowStyle: CSSProperties = { display: "flex", alignItems: "center", gap: 2 };
const rowNumberStyle: CSSProperties = { width: 32, textAlign: "right", paddingRight: 6, color: "#8796a8", fontSize: "0.72rem" };
const residueStyle = (picked: boolean, modified: boolean, terminus: boolean): CSSProperties => ({
  width: terminus ? 28 : 22,
  height: 22,
  marginLeft: terminus ? 38 : 0,
  padding: 0,
  fontFamily: "inherit",
  fontSize: "0.8rem",
  cursor: "pointer",
  borderRadius: 3,
  border: picked ? `2px solid ${SIM_COLOR}` : "1px solid #d5deeb",
  background: modified ? SIM_SELECTED_COLOR : terminus ? "#eef2f7" : "#fbfcfe",
  color: modified ? "#ffffff" : "#24364d",
  fontWeight: modified ? 700 : 400
});
const pickRowStyle: CSSProperties = { display: "flex", alignItems: "center", gap: 8, margin: "6px 0" };
const selectStyle: CSSProperties = { fontSize: "0.8rem", padding: "2px 4px", border: "1px solid #cbd6e5", borderRadius: 4 };
const modListStyle: CSSProperties = { listStyle: "none", padding: 0, margin: "4px 0" };
const modItemStyle: CSSProperties = {
  display: "flex",
  justifyContent: "space-between",
  alignItems: "center",
  padding: "2px 6px",
  borderBottom: "1px solid #eef2f7"
};
const massStyle: CSSProperties = { marginTop: 10, fontSize: "0.9rem" };
const pairStyle: CSSProperties = { display: "grid", gridTemplateColumns: "1fr 1fr", gap: 12 };
const cancelStyle: CSSProperties = {
  border: "1px solid #cbd6e5",
  background: "#f4f7fb",
  borderRadius: 4,
  cursor: "pointer",
  padding: "6px 12px",
  fontSize: "0.82rem"
};
const addStyle: CSSProperties = {
  border: `1px solid ${SIM_COLOR}`,
  background: SIM_COLOR,
  color: "#ffffff",
  borderRadius: 4,
  cursor: "pointer",
  padding: "6px 12px",
  fontWeight: 600,
  fontSize: "0.82rem"
};
const addDisabledStyle: CSSProperties = {
  background: "#e9edf3",
  borderColor: "#cbd6e5",
  color: "#7a8899",
  cursor: "default"
};
