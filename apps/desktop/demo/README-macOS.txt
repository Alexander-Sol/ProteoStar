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
The app opens with one minute of real Jurkat top-down data (rep2 fraction 7, 44-45 min,
22 MS1 scans), its 157 confident proteoform IDs (q <= 0.01), and a simulation fitted to
them (131 proteoforms).

- Left: an MS1 spectrum from the middle of the window. The real scan is drawn upward in
  blue, and the simulated scan of the same time is mirrored below in red.
- Right: the Simulation panel.

Things to try
-------------
- Drag "Noise level" and "Noise density" at the top of the panel. Untick Noise > enabled
  to see the clean simulation. The rest of the noise model is under "Advanced options".
- Click a proteoform in the panel's list. Its own envelope is drawn in amber over the
  simulated scan. Change its abundance, mass, RT or charge distribution; the simulated
  spectrum updates as you drag a slider.
- Change "Peak width k" under Run.
- Press the left/right arrow keys to move between scans. The simulated scan follows.
- Zoom with the zoom tool in the plot's top-right toolbar, then drag. Double-click to reset.
- "Save..." writes the edited simulation to a .simspec.json file.
