import type { CSSProperties, ReactNode } from "react";

export function WorkspaceBadge({ label }: { label: ReactNode }) {
  return <span style={badgeStyle}>{label}</span>;
}

export function StatusBanner({
  tone,
  children
}: {
  tone: "info" | "error" | "muted";
  children: ReactNode;
}) {
  return <div style={statusStyles[tone]}>{children}</div>;
}

export function MetricReadout({
  label,
  value
}: {
  label: string;
  value: ReactNode;
}) {
  return (
    <div style={metricStyle}>
      <span style={metricLabelStyle}>{label}</span>
      <span style={metricValueStyle}>{value}</span>
    </div>
  );
}

const badgeStyle: CSSProperties = {
  display: "inline-flex",
  alignItems: "center",
  padding: "0.5rem 0.75rem",
  borderRadius: "999px",
  border: "1px solid #c9d5e6",
  backgroundColor: "#f7faff",
  color: "#28415f",
  fontSize: "0.9rem",
  fontWeight: 600
};

const metricStyle: CSSProperties = {
  display: "grid",
  gap: 4,
  minWidth: 120
};

const metricLabelStyle: CSSProperties = {
  fontSize: "0.73rem",
  textTransform: "uppercase",
  letterSpacing: "0.08em",
  color: "#6a7d96"
};

const metricValueStyle: CSSProperties = {
  fontSize: "0.96rem",
  fontWeight: 600,
  color: "#203148"
};

const baseStatusStyle: CSSProperties = {
  borderRadius: 16,
  padding: "14px 16px",
  border: "1px solid #d8e2ef",
  fontSize: "0.95rem",
  lineHeight: 1.5
};

const statusStyles: Record<"info" | "error" | "muted", CSSProperties> = {
  info: {
    ...baseStatusStyle,
    backgroundColor: "#eef6ff",
    color: "#1d4067"
  },
  error: {
    ...baseStatusStyle,
    backgroundColor: "#fff1f2",
    border: "1px solid #ffcdd4",
    color: "#8a2131"
  },
  muted: {
    ...baseStatusStyle,
    backgroundColor: "#f7faff",
    color: "#506279"
  }
};
