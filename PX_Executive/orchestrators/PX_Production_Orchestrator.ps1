<#
    PREDATOR X :: LIVE PRODUCTION ORCHESTRATOR
    ARCHITECT: JAMES A. TILLAR
    STATUS:    ACTIVE // NO SIMULATION
    WINDOW:    3 HOURS (HARD STOP)
#>

# --- 1. CONFIGURATION: THE PHYSICS ---
$RunDurationHours = 3
$RootPath         = "E:\foundation"
$PythonExe        = "python" # Assumes active environment

# ASSETS (Validating against your verified repos)
$HarnessScript    = "$RootPath\mount-olympus-overnight-run\main_harness.py"
$MonitorScript    = "$RootPath\PX_Audit\Drift_Monitor.py"

# PIPELINE STAGING (The Staging Area for raw output)
$StagingArea      = "$RootPath\__PIPELINE_STAGING__"
$SurvivorsVault   = "$RootPath\SURVIVORS_Verified"
$GraveyardVault   = "$RootPath\GRAVEYARD_Failures"

# --- 2. INITIALIZATION: SECURING THE PERIMETER ---
Clear-Host
$StartTime = Get-Date
$EndTime   = $StartTime.AddHours($RunDurationHours)

# Set Console Physics
$Host.UI.RawUI.WindowTitle = "PREDATOR X :: PRODUCTION :: T-MINUS 03:00:00"

Write-Host "===================================================" -ForegroundColor Cyan
Write-Host "   PREDATOR X :: LIVE PRODUCTION EVENT" -ForegroundColor White
Write-Host "   ARCHITECT: JAMES TILLAR" -ForegroundColor Yellow
Write-Host "===================================================" -ForegroundColor Cyan
Write-Host "   START:  $($StartTime.ToString('HH:mm:ss'))" -ForegroundColor Gray
Write-Host "   CUTOFF: $($EndTime.ToString('HH:mm:ss'))" -ForegroundColor Red
Write-Host "===================================================" -ForegroundColor Cyan

# Create/Verify Physical Directories
New-Item -ItemType Directory -Force -Path $StagingArea | Out-Null
New-Item -ItemType Directory -Force -Path $SurvivorsVault | Out-Null
New-Item -ItemType Directory -Force -Path $GraveyardVault | Out-Null

# --- 3. SENTINEL ACTIVATION (Drift Monitor) ---
Write-Host ">>> [SENTINEL] Detaching Drift Monitor..." -ForegroundColor Gray
if (Test-Path $MonitorScript) {
    # Spawns separate process to watch the heartbeat independently
    Start-Process powershell -ArgumentList "-NoExit", "-Command", "& { $PythonExe '$MonitorScript'; Read-Host 'Monitor Halted' }"
    Write-Host "    + MONITOR ACTIVE" -ForegroundColor Green
} else {
    Write-Host "    ! CRITICAL: Monitor Script Missing. Proceeding Blind." -ForegroundColor Red
}

# --- 4. THE EXECUTION LOOP (The 3-Hour Raindrop) ---
Write-Host ">>> [PIPELINE] Engaging 35D Generation Sequence..." -ForegroundColor Yellow
$BatchCount = 0
$MoleculesProcessed = 0

while ((Get-Date) -lt $EndTime) {
    $Now = Get-Date
    $TimeRemaining = $EndTime - $Now
    $Host.UI.RawUI.WindowTitle = "PREDATOR X :: ACTIVE :: T-MINUS $($TimeRemaining.ToString('hh\:mm\:ss'))"

    # A. EXECUTE THE HARNESS (The Generator)
    # We call the python script to generate a batch into the Staging Area
    # Arguments: --mode production --target AUTO (derived from Dept 03) --out [Staging]
    Write-Host "[$($Now.ToString('HH:mm:ss'))] EXEC :: BATCH #$($BatchCount+1)" -NoNewline

    $ProcessInfo = New-Object System.Diagnostics.ProcessStartInfo
    $ProcessInfo.FileName = $PythonExe
    $ProcessInfo.Arguments = "`"$HarnessScript`" --mode production --output `"$StagingArea`""
    $ProcessInfo.RedirectStandardOutput = $true
    $ProcessInfo.UseShellExecute = $false
    $Process = [System.Diagnostics.Process]::Start($ProcessInfo)
    $Process.WaitForExit()

    # B. PROCESS THE OUTPUT (The Router)
    # Move files from Staging to their final destination based on the filename or content tags
    $RawFiles = Get-ChildItem -Path $StagingArea -Filter "*.json"

    if ($RawFiles.Count -gt 0) {
        foreach ($File in $RawFiles) {
            # READ THE VERDICT (Fail-Closed Logic)
            # We scan the first few bytes for the "status" field to avoid full parsing overhead
            $Content = Get-Content $File.FullName -Raw
            
            if ($Content -match '"status":\s*"SURVIVOR"' -or $Content -match '"verdict":\s*"PASS"') {
                Move-Item -Path $File.FullName -Destination $SurvivorsVault -Force
                $StatusColor = "Green"
                $Symbol = "+"
            } else {
                # Default to Graveyard (Fail-Closed)
                Move-Item -Path $File.FullName -Destination $GraveyardVault -Force
                $StatusColor = "DarkGray"
                $Symbol = "x"
            }
            $MoleculesProcessed++
        }
        Write-Host " [OK] ($($RawFiles.Count) Processed)" -ForegroundColor Green
    } else {
        Write-Host " [NO YIELD] (Resonance Tuning...)" -ForegroundColor DarkGray
    }

    $BatchCount++
    
    # NO SLEEP. The loop restarts immediately as requested.
}

# --- 5. MISSION DEBRIEF (End of Line) ---
$FinalTime = Get-Date
Write-Host "`n===================================================" -ForegroundColor Cyan
Write-Host "   PRODUCTION EVENT COMPLETE" -ForegroundColor Yellow
Write-Host "   TOTAL BATCHES: $BatchCount" -ForegroundColor Gray
Write-Host "   MOLECULES:     $MoleculesProcessed" -ForegroundColor Gray
Write-Host "   SURVIVORS:     $SurvivorsVault" -ForegroundColor Green
Write-Host "===================================================" -ForegroundColor Cyan