# Quick test orchestrator - 30 second run
$RunDurationHours = 0.0083  # 30 seconds
$GracePeriodSec   = 5   
$MaxFailures      = 3    

$RootPath         = "E:\foundation"
$PythonExe        = "python" 
$HarnessScript    = "$RootPath\PX_Audit\Autonomous_Research_Cycle.py"
$StagingArea      = "$RootPath\__PIPELINE_STAGING__"
$SurvivorsVault   = "$RootPath\SURVIVORS_Verified"
$LogVault         = "$RootPath\PX_LOGS"

# Create directories
$Directories = @($StagingArea, $SurvivorsVault, $LogVault)
foreach ($Dir in $Directories) { 
    if (-not (Test-Path $Dir)) { 
        New-Item -ItemType Directory -Force -Path $Dir | Out-Null 
    } 
}

Write-Host "[INFO] QUICK TEST: 30-second autonomous cycle" -ForegroundColor Cyan
Write-Host "[INFO] Harness: $HarnessScript" -ForegroundColor Gray

$StartTime = Get-Date
$EndTime   = $StartTime.AddHours($RunDurationHours)
$BatchCount = 0
$ConsecutiveFailures = 0

while (((Get-Date) -lt $EndTime) -and ($ConsecutiveFailures -lt $MaxFailures)) {
    Write-Host "`n[EXEC] Batch #$($BatchCount+1) starting..." -ForegroundColor Yellow
    
    $BatchLogOut = "$LogVault\TestBatch_$($BatchCount)_stdout.log"
    $BatchLogErr = "$LogVault\TestBatch_$($BatchCount)_stderr.log"

    # Fixed command line
    $CmdLine = "/c SET PYTHONPATH=E:\foundation && $PythonExe `"$HarnessScript`" > `"$BatchLogOut`" 2> `"$BatchLogErr`""
    
    $ProcInfo = New-Object System.Diagnostics.ProcessStartInfo
    $ProcInfo.FileName = "cmd.exe"
    $ProcInfo.Arguments = $CmdLine
    $ProcInfo.UseShellExecute = $false
    $ProcInfo.CreateNoWindow = $true

    $Harness = [System.Diagnostics.Process]::Start($ProcInfo)
    
    if (-not $Harness.WaitForExit(30000)) {
        Write-Host "[ERROR] Timeout!" -ForegroundColor Red
        $Harness.Kill()
        $ConsecutiveFailures++
    } else {
        if ($Harness.ExitCode -ne 0) {
            $ErrContent = Get-Content $BatchLogErr -Raw -ErrorAction SilentlyContinue
            Write-Host "[ERROR] Harness crash (Exit: $($Harness.ExitCode))" -ForegroundColor Red
            Write-Host "  Error: $ErrContent" -ForegroundColor DarkRed
            $ConsecutiveFailures++
        } else {
            Write-Host "[SUCCESS] Batch completed!" -ForegroundColor Green
            $ConsecutiveFailures = 0
            
            # Check if worldline was created
            $NewWorldline = Get-ChildItem "$RootPath\PX_Warehouse\WorldLines" | 
                Where-Object { $_.LastWriteTime -gt $StartTime } | 
                Select-Object -Last 1
            
            if ($NewWorldline) {
                Write-Host "  Created: $($NewWorldline.Name)" -ForegroundColor Green
            }
        }
    }
    
    $BatchCount++
    Start-Sleep -Seconds 2
}

Write-Host "`n========================================" -ForegroundColor Cyan
Write-Host "TEST COMPLETE" -ForegroundColor Yellow
Write-Host "  Batches executed: $BatchCount" -ForegroundColor White
Write-Host "  Consecutive failures: $ConsecutiveFailures" -ForegroundColor White
Write-Host "========================================" -ForegroundColor Cyan
