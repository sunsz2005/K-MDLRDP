$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $PSScriptRoot
$instances = Join-Path $projectRoot "Instances\PPI"
$output = Join-Path $projectRoot "Results\formal\ppi_gate_greedy.csv"
Set-Location $projectRoot
& "python" .\python\run_experiments.py $instances `
    --methods greedy --k "1,2,3,4,5" --seeds 13579 `
    --time-limit 3600 --output $output
