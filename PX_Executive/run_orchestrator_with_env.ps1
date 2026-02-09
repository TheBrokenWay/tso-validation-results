# PREDATOR X — Industrial loop with correct Conda environment
# Run this in PowerShell after stopping any background orchestrator (Ctrl+C).
# Requires: numpy and rdkit in the .conda env (pip install numpy rdkit).

$ErrorActionPreference = "Stop"
Set-Location E:\foundation

# Use conda run so the env is guaranteed (no activate persistence issues)
$env:PRV_LIVE_RESEARCH = "1"
conda run -p E:\foundation\.conda python PX_Executive/PRV_24H_Orchestrator.py *> PX_LOGS/live_run_output.l
