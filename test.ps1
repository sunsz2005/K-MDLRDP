param([string]$Python = "python")

$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
& (Join-Path $projectRoot "build.ps1") -Configuration Debug
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

& (Join-Path $projectRoot "bin\test_core.exe")
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

& $Python (Join-Path $projectRoot "tests\test_exact.py")
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

& $Python (Join-Path $projectRoot "tests\test_hybrid.py")
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host "All tests passed"

