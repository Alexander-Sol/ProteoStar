// The simulation workbench's shared look: its colours and the input widgets its panel and
// dialogs are built from.

import { useEffect, useRef, useState } from "react";
import type { CSSProperties, ReactNode } from "react";

/** The simulated run's colour wherever it is drawn (mirror traces, panel accents). */
export const SIM_COLOR = "#b5482f";
/** The selected proteoform's own envelope, drawn over the simulated scan. */
export const SIM_SELECTED_COLOR = "#e8a317";

export function SmallButton({
  onClick,
  title,
  children
}: {
  onClick: () => void;
  title?: string;
  children: ReactNode;
}) {
  return (
    <button type="button" onClick={onClick} title={title} style={smallButtonStyle}>
      {children}
    </button>
  );
}

export function TextField({
  label,
  value,
  onChange
}: {
  label: string;
  value: string;
  onChange: (v: string) => void;
}) {
  return (
    <label style={fieldStyle}>
      <span style={fieldLabelStyle}>{label}</span>
      <input style={inputStyle} value={value} onChange={(e) => onChange(e.target.value)} />
    </label>
  );
}

/**
 * A number input, optionally with a slider (linear or log). The text is a local draft: a value
 * that does not parse to a finite number is shown but never sent, so a half-typed "1e" does not
 * reach the simulator. Values that parse but the backend rejects (σ ≤ 0, min > max) come back as
 * the panel's error line.
 */
export function NumberField({
  label,
  hint,
  value,
  step,
  decimals = 2,
  slider,
  onChange
}: {
  label: string;
  hint?: string;
  value: number;
  step?: number;
  /** Decimal places shown in the box; the stored value keeps full precision. */
  decimals?: number;
  slider?: { min: number; max: number; log?: boolean };
  onChange: (v: number) => void;
}) {
  const format = (v: number): string =>
    // Values too small for the fixed places (peak width k ~1e-7) keep three significant figures.
    v !== 0 && Math.abs(v) < 10 ** -decimals ? v.toPrecision(3) : String(Number(v.toFixed(decimals)));
  const [draft, setDraft] = useState(format(value));
  const focused = useRef(false);
  useEffect(() => {
    // Follow outside changes (slider, undo-by-load) unless the user is mid-edit here.
    if (!focused.current) setDraft(format(value));
  }, [value]);

  const SLIDER_STEPS = 1000;
  const toSlider = (v: number): number => {
    if (!slider) return 0;
    const { min, max, log } = slider;
    const f = log
      ? (Math.log(Math.max(v, min)) - Math.log(min)) / (Math.log(max) - Math.log(min))
      : (v - min) / (max - min);
    return Math.round(Math.max(0, Math.min(1, f)) * SLIDER_STEPS);
  };
  const fromSlider = (pos: number): number => {
    const { min, max, log } = slider!;
    const f = pos / SLIDER_STEPS;
    const v = log ? Math.exp(Math.log(min) + f * (Math.log(max) - Math.log(min))) : min + f * (max - min);
    // Snap to the field's step; otherwise three significant figures, since slider positions are
    // approximate by nature.
    return step ? Number((Math.round(v / step) * step).toFixed(10)) : Number(v.toPrecision(3));
  };

  return (
    <label style={fieldStyle}>
      <span style={fieldLabelStyle}>
        {label}
        {hint ? <span style={hintStyle}> {hint}</span> : null}
      </span>
      <input
        style={inputStyle}
        type="number"
        step={step ?? "any"}
        value={draft}
        onFocus={() => (focused.current = true)}
        onBlur={() => {
          focused.current = false;
          setDraft(format(value));
        }}
        onChange={(e) => {
          setDraft(e.target.value);
          const v = Number(e.target.value);
          if (e.target.value.trim() !== "" && Number.isFinite(v)) onChange(v);
        }}
      />
      {slider ? (
        <input
          type="range"
          min={0}
          max={SLIDER_STEPS}
          value={toSlider(value)}
          onChange={(e) => onChange(fromSlider(Number(e.target.value)))}
          style={sliderStyle}
        />
      ) : null}
    </label>
  );
}

const smallButtonStyle: CSSProperties = {
  border: "1px solid #cbd6e5",
  background: "#f4f7fb",
  borderRadius: 4,
  cursor: "pointer",
  padding: "2px 8px",
  fontSize: "0.76rem"
};
const fieldStyle: CSSProperties = { display: "flex", flexDirection: "column", gap: 2, margin: "4px 0" };
const fieldLabelStyle: CSSProperties = { color: "#3a4a60" };
const hintStyle: CSSProperties = { color: "#8796a8" };
const inputStyle: CSSProperties = {
  fontSize: "0.78rem",
  padding: "2px 4px",
  border: "1px solid #cbd6e5",
  borderRadius: 4,
  width: "100%",
  boxSizing: "border-box"
};
const sliderStyle: CSSProperties = { width: "100%", accentColor: SIM_COLOR };
