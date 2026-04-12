import type { ChangeEvent } from "react";

export function FileOpenButton({
  onSelect
}: {
  onSelect(file: File): void;
}) {
  return (
    <label style={labelStyle}>
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
      Open `.imsp` File
    </label>
  );
}

const labelStyle = {
  display: "inline-flex",
  alignItems: "center",
  justifyContent: "center",
  gap: 8,
  minHeight: 44,
  borderRadius: 999,
  padding: "0 1rem",
  border: "1px solid #2f5f95",
  background: "linear-gradient(180deg, #2f6fb0 0%, #265c92 100%)",
  color: "#ffffff",
  fontWeight: 700,
  cursor: "pointer",
  whiteSpace: "nowrap"
} satisfies Record<string, string | number>;
