$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $PSScriptRoot
$instances = Join-Path $projectRoot "Instances\PPI"
$output = Join-Path $projectRoot "Results\formal\ppi_formal_ga_g_hybrid_k1_extra.csv"

Set-Location $projectRoot
& "python" .\python\run_experiments.py $instances `
  --methods "ga-g-hybrid" `
  --k "1" `
  --seeds "24680,97531,86420,112358,271828,314159,161803,424242,987654" `
  --time-limit 3600 `
  --output $output

if ($LASTEXITCODE -ne 0) {
  exit $LASTEXITCODE
}
