MsViewer simulation workbench demo
==================================

Opening it
----------
1. Unzip and drag "MsViewer Demo.app" to Applications (or run it from anywhere).
2. The app is not notarized, so the first launch is blocked. Either:
   - right-click the app > Open > Open, or, if macOS still refuses,
     System Settings > Privacy & Security > "Open Anyway"; or
   - in Terminal:  xattr -cr "/Applications/MsViewer Demo.app"
   After the first launch it opens normally.

What loads
----------
The app opens with ten minutes of real Jurkat top-down data (rep2 fraction 7, 40-50 min,
232 MS1 scans), its 822 confident proteoform IDs (q <= 0.01), and a simulation fitted to
them (337 proteoforms).

- Left: an MS1 spectrum from the middle of the window. The real scan is drawn upward in
  blue, and the simulated scan of the same time is mirrored below in red.
  "Chromatogram" (top of the plot) switches to the whole run: the real TIC above the
  simulated TIC, mirrored, with the selected proteoform's XIC in amber in both.
  "Spectrum" switches back.
- Right: the Simulation panel.

Things to try
-------------
- Drag "Noise level" and "Noise density" at the top of the panel. Untick Noise > enabled
  to see the clean simulation. The rest of the noise model is under "Advanced options".
- Click a proteoform in the panel's list. Its own envelope is drawn in amber over the
  simulated scan. Change its abundance, retention time, RT σ (chromatographic peak width)
  or mean charge; the simulated spectrum updates as you drag a slider. The elution tail and
  charge spread are under "Advanced options".
- In the chromatogram, drag a proteoform's abundance or retention time and watch the
  simulated TIC and its XIC move. The real XIC is extracted at the simulated proteoform's
  strongest peaks, so it changes when its mean charge does. "Zoom to proteoform" scales the
  axis to the XICs. Click the chromatogram to pick the scan the spectrum view shows.
- "+ Add proteoform" opens a dialog: paste a sequence (or a FASTA entry), click residues to
  add modifications, and the monoisotopic mass is computed for you. "Predict charge from
  sequence" sets the mean charge from mass and composition (hover it for the model).
- Change "Peak width k" under Run.
- Press the left/right arrow keys to move between scans. The simulated scan follows.
- Zoom with the zoom tool in the plot's top-right toolbar, then drag. Double-click to reset.
- "Save..." writes the edited simulation to a .simspec.json file.
