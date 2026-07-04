<#
.SYNOPSIS
    Ralph loop: repeatedly spin up fresh headless Claude Code instances against a plan file.

.DESCRIPTION
    Each iteration launches a NEW `claude` process (fresh context window). The agent is told to:
      1. Read the plan file.
      2. Pick the single next incomplete task and do it.
      3. Write progress back into the plan (mark done / leave notes for the next instance).
      4. Stop after ONE task -- or, if context is filling up, save progress and stop early.

    Durable state lives in the PLAN FILE, not in the model's context. "Context window fills up"
    is handled structurally: every iteration is a brand-new instance, so each task gets a full
    fresh window; the agent only has to flush progress to the plan before it exits.

    The loop ends when the agent emits the sentinel ALL_TASKS_COMPLETE, when MaxIterations is
    hit, or after too many consecutive failures.

.EXAMPLE
    .\ralph-loop.ps1 -PlanPath .\agent_info\FlashLFQ-Rust-Rewrite-Feasibility.md
#>

[CmdletBinding()]
param(
    # Path to the plan / task-list markdown that the agent reads and updates each iteration.
    [string]$PlanPath = ".\PLAN.md",

    # Safety cap on iterations so a misbehaving loop can't run forever.
    [int]$MaxIterations = 25,

    # Bail out after this many consecutive failed `claude` invocations.
    [int]$MaxConsecutiveFailures = 3,

    # Seconds to wait between iterations.
    [int]$SleepSeconds = 5,

    # Optional model override (e.g. "opus", "sonnet"). Empty = account default.
    [string]$Model = "",

    # Directory to run in (repo root by default = the script's own folder).
    [string]$WorkDir = $PSScriptRoot
)

$ErrorActionPreference = "Stop"
Set-Location $WorkDir

if (-not (Test-Path $PlanPath)) {
    throw "Plan file not found: $PlanPath"
}
$resolvedPlan = (Resolve-Path $PlanPath).Path

# Per-run log directory.
$logDir = Join-Path $WorkDir ".ralph-logs"
if (-not (Test-Path $logDir)) { New-Item -ItemType Directory -Path $logDir | Out-Null }

# The instruction handed to each fresh instance. Single durable task per run.
$prompt = @"
You are one iteration of an automated loop. A plan file lives at: $resolvedPlan

Do exactly this:
1. Read the plan file in full.
2. Identify the SINGLE next incomplete task.
3. Implement that one task completely and correctly. Build/test as appropriate.
4. Update the plan file to record what you did: mark the task complete and leave a short
   'next up' note so the following instance has the context it needs.
5. Commit your work with a clear message (you are on a working branch; do not push).
6. Then STOP. Do not start a second task.

If your context window is filling up before the task is finished, do NOT push forward blindly:
write your partial progress and a precise resume note into the plan file, commit, and stop.

If EVERY task in the plan is already complete, make no changes and output the single token:
ALL_TASKS_COMPLETE
"@

# Build the claude argument list.
$claudeArgs = @("-p", $prompt, "--dangerously-skip-permissions")
if ($Model -ne "") { $claudeArgs += @("--model", $Model) }

$consecutiveFailures = 0

for ($i = 1; $i -le $MaxIterations; $i++) {
    $stamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
    Write-Host "=== Ralph iteration $i / $MaxIterations  ($stamp) ===" -ForegroundColor Cyan

    $iterLog = Join-Path $logDir ("iter-{0:D3}.log" -f $i)

    # Fresh instance every iteration. Capture stdout so we can scan for the sentinel.
    $output = & claude @claudeArgs 2>&1 | Tee-Object -FilePath $iterLog
    $exit = $LASTEXITCODE

    $outText = ($output | Out-String)

    if ($exit -ne 0) {
        $consecutiveFailures++
        Write-Host "claude exited with code $exit (failure $consecutiveFailures/$MaxConsecutiveFailures). See $iterLog" -ForegroundColor Yellow
        if ($consecutiveFailures -ge $MaxConsecutiveFailures) {
            Write-Host "Too many consecutive failures. Stopping." -ForegroundColor Red
            break
        }
        Start-Sleep -Seconds $SleepSeconds
        continue
    }
    $consecutiveFailures = 0

    if ($outText -match "ALL_TASKS_COMPLETE") {
        Write-Host "Agent reports all tasks complete. Loop finished after $i iteration(s)." -ForegroundColor Green
        break
    }

    Write-Host "Iteration $i done. Pausing $SleepSeconds s before next instance..." -ForegroundColor DarkGray
    Start-Sleep -Seconds $SleepSeconds
}

Write-Host "Ralph loop exited. Logs in: $logDir" -ForegroundColor Cyan
