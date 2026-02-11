# Nipah Analysis — entry script (validation gate then summary)
# Control flow: Load manifest → Validate manifest → Discover raw → Validate each file (checksum, schema, strain, CFR, temporal) → FAIL if any reject → Proceed

$ErrorActionPreference = "Stop"
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $ScriptDir

$RawDir   = Join-Path $ScriptDir "data\raw"
$NormDir  = Join-Path $ScriptDir "data\normalized"
$ResultsDir = Join-Path $ScriptDir "results"

foreach ($d in $RawDir, $NormDir, $ResultsDir) {
    if (-not (Test-Path $d)) { New-Item -ItemType Directory -Path $d -Force | Out-Null }
}

# Blocking validation gate: Python validator (manifest + raw data)
# Pass --allow-multi-strain if raw data contains more than one strain (NiV_Malaysia + NiV_Bangladesh)
$pythonScript = Join-Path $ScriptDir "run_nipah_analysis.py"
if (-not (Test-Path $pythonScript)) {
    Write-Error "run_nipah_analysis.py not found at $pythonScript"
    exit 1
}
$args = $args  # forward any args (e.g. --allow-multi-strain)
& python $pythonScript @args
if ($LASTEXITCODE -ne 0) {
    Write-Error "Nipah validation FAILED. No partial runs. Fix errors and re-run."
    exit $LASTEXITCODE
}

Write-Host "Nipah analysis complete. See results\summary.json and results\validation_report.json"
