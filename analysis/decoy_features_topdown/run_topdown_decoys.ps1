# Top-down decoy feature-detection experiment.
# Detect target + 3 decoys on a top-down raw, then score each with the 2-D/3-D envelope-fit cosine.
# All runs share one config so the envelope model is the ONLY difference:
#   TOPDOWN=1 APEX_PREGATE=0 REFINE_METHOD=shift_apex TD_MONO_FIT=0 ISODEC_CHARGE=0
# shift_apex refine is decoy-aware (Deconvoluter bound to the run's EnvelopeModel); the multi refine,
# averagine TD_MONO_FIT, and IsoDec are disabled because they would re-anchor a decoy onto real
# averagine patterns (laundering). Detection uses the decoy-threaded multi-charge ladder.
#
# Usage: run_topdown_decoys.ps1 <dataset-key> <raw-path>
param(
  [Parameter(Mandatory=$true)][string]$Ds,
  [Parameter(Mandatory=$true)][string]$Raw
)
$ErrorActionPreference = "Continue"   # native exe stderr must NOT abort the run (PS 5.1 NativeCommandError)
$here = "F:\ProteoStar\analysis\decoy_features_topdown"
$det  = "F:\ProteoStar\target\release\examples\detect_features_tsv.exe"
$scr  = "F:\ProteoStar\target\release\examples\decoy_score_export.exe"

$models = @(
  @{ name="target";   env=@{};                                                          smodel="averagine"; scale="1.0"    },
  @{ name="shifted";  env=@{ DECOY_LATTICE="scaled"; DECOY_SPACING_SCALE="0.9368" };     smodel="averagine"; scale="0.9368" },
  @{ name="weird";    env=@{ COMB_MODEL="custom" };                                      smodel="custom";    scale="1.0"    },
  @{ name="shuffled"; env=@{ COMB_MODEL="shuffled" };                                    smodel="shuffled";  scale="1.0"    }
)
$base = @{ TOPDOWN="1"; APEX_PREGATE="0"; REFINE_METHOD="shift_apex"; TD_MONO_FIT="0"; ISODEC_CHARGE="0"; DETECT_PROGRESS="1" }
$decoyKeys = @("COMB_MODEL","DECOY_LATTICE","DECOY_SPACING_SCALE")

foreach ($m in $models) {
  $name = $m.name
  $out  = Join-Path $here "t_${Ds}_${name}.tsv"
  $dlog = Join-Path $here "log_${Ds}_${name}_detect.txt"
  Write-Output "==== [$Ds/$name] DETECT $(Get-Date -Format HH:mm:ss) ===="
  foreach ($k in $decoyKeys) { Remove-Item "env:$k" -ErrorAction SilentlyContinue }
  foreach ($k in $base.Keys) { Set-Item "env:$k" $base[$k] }
  foreach ($k in $m.env.Keys) { Set-Item "env:$k" $m.env[$k] }
  $sw = [System.Diagnostics.Stopwatch]::StartNew()
  & $det $Raw $out *> $dlog
  $sw.Stop()
  Write-Output "  detect wall: $([math]::Round($sw.Elapsed.TotalSeconds,1))s  (log: $dlog)"
  Get-Content $dlog | Select-Object -Last 3 | ForEach-Object { Write-Output "    $_" }

  $refined = Join-Path $here "t_${Ds}_${name}.refined.tsv"
  $sout    = Join-Path $here "s_${Ds}_${name}.tsv"
  $slog    = Join-Path $here "log_${Ds}_${name}_score.txt"
  Write-Output "---- [$Ds/$name] SCORE ($($m.smodel) $($m.scale)) ----"
  $env:SCORE_MAX_ISOTOPES = "60"
  $env:SCORE_RT_SIGMA     = "0.30"
  & $scr $Raw $refined $sout $m.smodel $m.scale *> $slog
  Get-Content $slog | Select-Object -Last 2 | ForEach-Object { Write-Output "    $_" }
}
Write-Output "==== [$Ds] DONE $(Get-Date -Format HH:mm:ss) ===="
