$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $PSScriptRoot
$instances = Join-Path $projectRoot "Instances\PPI"
$output = Join-Path $projectRoot "Results\formal\ppi_5layer_k5_vns_extra.csv"

Set-Location $projectRoot
& "python" .\python\run_experiments.py $instances `
  --include-regex "PPI-5Layers$" `
  --methods "vns-g-fi" `
  --k "5" `
  --seeds "24680,97531,86420,112358,271828,314159,161803,424242,987654" `
  --time-limit 3600 `
  --output $output

if ($LASTEXITCODE -ne 0) {
  exit $LASTEXITCODE
}
