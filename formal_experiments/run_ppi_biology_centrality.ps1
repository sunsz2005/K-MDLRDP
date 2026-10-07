$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $PSScriptRoot
$biology = Join-Path $projectRoot "Results\biology"
$instance = Join-Path $projectRoot "Instances\PPI\PPI-5Layers"
$executable = Join-Path $projectRoot "bin\ppi_exact_centrality.exe"

New-Item -ItemType Directory -Path $biology -Force | Out-Null
Set-Location $projectRoot
& "python" .\python\prepare_ppi_biology.py
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

for ($layer = 1; $layer -le 5; $layer++) {
  Write-Output "RUN exact centrality layer $layer"
  & $executable (Join-Path $instance "layer_$layer.txt") `
    (Join-Path $biology "ppi5_layer_${layer}_centrality.csv")
  if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
}
