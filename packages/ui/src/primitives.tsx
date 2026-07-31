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

/** An indeterminate loading spinner. Self-contained SVG using SMIL `animateTransform`, so it needs
 *  no global `@keyframes` / CSS file and animates even inside the system webview. */
export function Spinner({
  size = 30,
  color = "#3b6ea5"
}: {
  size?: number;
  color?: string;
}) {
  const c = size / 2;
  const r = c - 3;
  const circumference = 2 * Math.PI * r;
  return (
    <svg
      width={size}
      height={size}
      viewBox={`0 0 ${size} ${size}`}
      role="status"
      aria-label="Loading"
    >
      <circle cx={c} cy={c} r={r} fill="none" stroke="#d8e2ef" strokeWidth={3} />
      <circle
        cx={c}
        cy={c}
        r={r}
        fill="none"
        stroke={color}
        strokeWidth={3}
        strokeLinecap="round"
        strokeDasharray={`${circumference * 0.28} ${circumference}`}
      >
        <animateTransform
          attributeName="transform"
          type="rotate"
          from={`0 ${c} ${c}`}
          to={`360 ${c} ${c}`}
          dur="0.9s"
          repeatCount="indefinite"
        />
      </circle>
    </svg>
  );
}

/** A centered spinner with an optional caption, sized to fill its container — the empty-state to
 *  drop into a plot/panel body while work is in flight. */
export function LoadingState({ message }: { message?: ReactNode }) {
  return (
    <div style={loadingStateStyle}>
      <Spinner />
      {message ? <span style={loadingMessageStyle}>{message}</span> : null}
    </div>
  );
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

const loadingStateStyle: CSSProperties = {
  display: "flex",
  flexDirection: "column",
  alignItems: "center",
  justifyContent: "center",
  gap: 12,
  width: "100%",
  height: "100%",
  minHeight: 120,
  color: "#506279",
  fontSize: "0.95rem",
  textAlign: "center"
};

const loadingMessageStyle: CSSProperties = {
  maxWidth: 340,
  lineHeight: 1.5
};

const metricStyle: CSSProperties = {
  display: "grid",
  gap: 2,
  minWidth: 80
};

const metricLabelStyle: CSSProperties = {
  fontSize: "0.62rem",
  textTransform: "uppercase",
  letterSpacing: "0.08em",
  color: "#6a7d96"
};

const metricValueStyle: CSSProperties = {
  fontSize: "0.78rem",
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
