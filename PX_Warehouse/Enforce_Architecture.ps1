# PILOT: ENFORCE_ARCHITECTURE.ps1
# PURPOSE: Destroys unauthorized folders in CommercialAssets
$RootPath = 'E:\foundation\PX_Warehouse\CommercialAssets'
$Mandatory = @('Audit_Trail', 'Dossier_Final', 'Executive_Summary', 'Gold', 'Silver', 'Learning_Material')

if (!(Test-Path $RootPath)) { New-Item -ItemType Directory -Path $RootPath | Out-Null }
foreach ($Folder in $Mandatory) {
    $Target = "$RootPath\$Folder"
    if (!(Test-Path $Target)) { New-Item -ItemType Directory -Force -Path $Target | Out-Null; Write-Host "Restored: $Folder" -ForegroundColor Green }
}

Get-ChildItem -Path $RootPath -Directory | ForEach-Object {
    if ($Mandatory -notcontains $_.Name) { Remove-Item $_.FullName -Recurse -Force; Write-Host "Purged: $($_.Name)" -ForegroundColor Red }
}
