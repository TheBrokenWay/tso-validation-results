# 🔄 The "Zombie" Loop: Revives the Miner every time it dies
# Run from repo root, or this script will cd to the repo root (parent of PX_Executive).

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$RepoRoot = Split-Path -Parent $ScriptDir
Set-Location $RepoRoot

# Clear any per-run limits so the miner processes the full queue each cycle
$env:PRV_MAX_ITEMS = ""

while ($true) {
    Write-Host "`n⛏️  MINER STARTING..." -ForegroundColor Cyan
    python PX_Executive/run_prv_repurposed.py
    Write-Host "`n😴 Resting for 5 seconds before restart..." -ForegroundColor Yellow
    Start-Sleep -Seconds 5
}
