import type { ChangeEvent, CSSProperties } from "react";

export function FileOpenButton({
  onSelect,
  label = "Open `.imsp` File",
  color = "#2f6fb0"
}: {
  onSelect(file: File): void;
  label?: string;
  color?: string;
}) {
  const style: CSSProperties = {
    display: "inline-flex",
    alignItems: "center",
    justifyContent: "center",
    gap: 8,
    minHeight: 44,
    borderRadius: 999,
    padding: "0 1rem",
    border: `1px solid ${color}`,
    background: `linear-gradient(180deg, ${color} 0%, ${darken(color)} 100%)`,
    color: "#ffffff",
    fontWeight: 700,
    cursor: "pointer",
    whiteSpace: "nowrap" as const
  };

  return (
    <label style={style}>
      <input
        accept=".imsp"
        onChange={(event: ChangeEvent<HTMLInputElement>) => {
          const file = event.target.files?.[0];
          if (file) {
            onSelect(file);
            event.target.value = "";
          }
        }}
        style={{ display: "none" }}
        type="file"
      />
      {label}
    </label>
  );
}

/** Darken a hex color by ~15% for the gradient bottom stop. */
function darken(hex: string): string {
  const n = parseInt(hex.replace("#", ""), 16);
  const r = Math.max(0, ((n >> 16) & 0xff) - 38);
  const g = Math.max(0, ((n >> 8) & 0xff) - 38);
  const b = Math.max(0, (n & 0xff) - 38);
  return `#${[r, g, b].map((c) => c.toString(16).padStart(2, "0")).join("")}`;
}
