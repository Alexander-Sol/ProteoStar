import { IMSP_MAGIC } from "@msbrowser/imsp-core";
import { PLOT_ADAPTER_NAME } from "@msbrowser/plot-adapter";
import { WorkspaceBadge } from "@msbrowser/ui";
import { createInitialViewerState } from "@msbrowser/viewer-state";

const viewerState = createInitialViewerState();

export default function HomePage() {
  return (
    <main className="page-shell">
      <section className="hero-card">
        <p className="eyebrow">MsBrowser</p>
        <h1>Workspace scaffold is ready for IMSP parsing and viewer work.</h1>
        <p className="lead">
          The repository now has isolated packages for parsing, viewer state,
          plotting, and UI so the web app can stay decoupled from binary file
          details.
        </p>

        <div className="badge-row">
          <WorkspaceBadge label={`magic ${IMSP_MAGIC}`} />
          <WorkspaceBadge label={viewerState.activePanel} />
          <WorkspaceBadge label={PLOT_ADAPTER_NAME} />
        </div>
      </section>
    </main>
  );
}
