import type { CSSProperties, ReactNode } from "react";

export function ViewerShell({
  title,
  subtitle,
  toolbar,
  children
}: {
  title: string;
  subtitle: string;
  toolbar?: ReactNode;
  children: ReactNode;
}) {
  return (
    <main style={shellStyle}>
      <header style={headerStyle}>
        <div>
          <p style={eyebrowStyle}>MsBrowser</p>
          <h1 style={titleStyle}>{title}</h1>
          <p style={subtitleStyle}>{subtitle}</p>
        </div>
        {toolbar ? <div style={toolbarStyle}>{toolbar}</div> : null}
      </header>
      <section style={contentStyle}>{children}</section>
    </main>
  );
}

export function Panel({
  header,
  children
}: {
  header: ReactNode;
  children: ReactNode;
}) {
  return (
    <section style={panelStyle}>
      {header}
      <div style={panelBodyStyle}>{children}</div>
    </section>
  );
}

export function PanelHeader({
  title,
  subtitle,
  actions,
  readouts
}: {
  title: string;
  subtitle?: string;
  actions?: ReactNode;
  readouts?: ReactNode;
}) {
  return (
    <div style={panelHeaderStyle}>
      <div style={panelTitleWrapStyle}>
        <div>
          <h2 style={panelTitleStyle}>{title}</h2>
          {subtitle ? <p style={panelSubtitleStyle}>{subtitle}</p> : null}
        </div>
        {readouts ? <div style={readoutRowStyle}>{readouts}</div> : null}
      </div>
      {actions ? <div style={actionsStyle}>{actions}</div> : null}
    </div>
  );
}

export function PanelActionButton({
  pressed = false,
  onClick,
  children
}: {
  pressed?: boolean;
  onClick?: () => void;
  children: ReactNode;
}) {
  return (
    <button onClick={onClick} style={buttonStyle(pressed)} type="button">
      {children}
    </button>
  );
}

const shellStyle: CSSProperties = {
  minHeight: "100vh",
  padding: "16px clamp(12px, 2vw, 24px) 20px",
  display: "grid",
  gap: 12
};

const headerStyle: CSSProperties = {
  display: "flex",
  alignItems: "flex-start",
  justifyContent: "space-between",
  gap: 16,
  flexWrap: "wrap"
};

const eyebrowStyle: CSSProperties = {
  margin: 0,
  fontSize: "0.76rem",
  textTransform: "uppercase",
  letterSpacing: "0.12em",
  color: "#4b6486",
  fontWeight: 700
};

const titleStyle: CSSProperties = {
  margin: "4px 0 0",
  fontSize: "clamp(1.5rem, 3vw, 2.1rem)",
  lineHeight: 1.05,
  color: "#13253d"
};

const subtitleStyle: CSSProperties = {
  margin: "8px 0 0",
  fontSize: "0.92rem",
  lineHeight: 1.45,
  maxWidth: 560,
  color: "#506279"
};

const toolbarStyle: CSSProperties = {
  display: "flex",
  alignItems: "center",
  gap: 12,
  flexWrap: "wrap"
};

const contentStyle: CSSProperties = {
  display: "grid",
  gap: 12
};

const panelStyle: CSSProperties = {
  display: "grid",
  gridTemplateRows: "auto minmax(380px, 1fr)",
  backgroundColor: "rgba(255,255,255,0.94)",
  border: "2px solid #d2ddea",
  borderRadius: 18,
  overflow: "hidden",
  boxShadow: "0 10px 28px rgba(22, 32, 51, 0.05)"
};

const panelBodyStyle: CSSProperties = {
  minHeight: 380,
  padding: 4
};

const panelHeaderStyle: CSSProperties = {
  display: "flex",
  justifyContent: "space-between",
  gap: 16,
  padding: "12px 16px 10px",
  borderBottom: "1px solid #e2e9f3",
  alignItems: "flex-start",
  flexWrap: "wrap"
};

const panelTitleWrapStyle: CSSProperties = {
  display: "grid",
  gap: 8
};

const panelTitleStyle: CSSProperties = {
  margin: 0,
  fontSize: "0.98rem",
  color: "#203148"
};

const panelSubtitleStyle: CSSProperties = {
  margin: "4px 0 0",
  color: "#61738c",
  fontSize: "0.85rem"
};

const readoutRowStyle: CSSProperties = {
  display: "flex",
  gap: 16,
  flexWrap: "wrap"
};

const actionsStyle: CSSProperties = {
  display: "flex",
  gap: 10,
  flexWrap: "wrap"
};

function buttonStyle(pressed: boolean): CSSProperties {
  return {
    appearance: "none",
    border: `1px solid ${pressed ? "#4c51bf" : "#cdd8e7"}`,
    backgroundColor: pressed ? "#eef0ff" : "#f9fbfe",
    color: pressed ? "#3730a3" : "#30445f",
    borderRadius: 999,
    padding: "0.55rem 0.9rem",
    fontSize: "0.9rem",
    fontWeight: 600,
    cursor: "pointer"
  };
}
