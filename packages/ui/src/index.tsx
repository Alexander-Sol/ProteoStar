import type { ReactNode } from "react";

export function WorkspaceBadge({ label }: { label: ReactNode }) {
  return <span style={badgeStyle}>{label}</span>;
}

const badgeStyle = {
  display: "inline-flex",
  alignItems: "center",
  padding: "0.5rem 0.75rem",
  borderRadius: "999px",
  border: "1px solid #c9d5e6",
  backgroundColor: "#f7faff",
  color: "#28415f",
  fontSize: "0.9rem",
  fontWeight: 600
} satisfies Record<string, string | number>;
