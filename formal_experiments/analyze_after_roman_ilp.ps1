param(
    [Parameter(Mandatory = $true)][int]$DirectProcessId,
    [Parameter(Mandatory = $true)][int]$CumulativeProcessId
)

$ErrorActionPreference = "Stop"
Wait-Process -Id $DirectProcessId -ErrorAction SilentlyContinue
Wait-Process -Id $CumulativeProcessId -ErrorAction SilentlyContinue

$root = Split-Path -Parent $PSScriptRoot
$python = "python"
Push-Location $root
try {
    & $python ".\python\prepare_roman_ilp_consensus.py"
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    & $python ".\python\analyze_roman_double_comparison.py"
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    & $python ".\python\analyze_degree_matched_go.py"
    exit $LASTEXITCODE
}
finally {
    Pop-Location
}
