$projectRoot = Split-Path -Parent $PSScriptRoot
$resultFile = Join-Path $projectRoot "Results\formal\ppi_5layer_k5_vns_extra.csv"
$pidFile = Join-Path $projectRoot "Results\formal\ppi_5layer_k5_vns.pid"

$rows = 0
if (Test-Path -LiteralPath $resultFile) {
  $rows = @(Import-Csv -LiteralPath $resultFile).Count
}

$pidValue = $null
$process = $null
if (Test-Path -LiteralPath $pidFile) {
  $pidValue = (Get-Content -LiteralPath $pidFile -Raw).Trim()
  if ($pidValue) {
    $process = Get-Process -Id ([int]$pidValue) -ErrorAction SilentlyContinue
  }
}

[pscustomobject]@{
  completed_new_runs = $rows
  expected_new_runs = 9
  percent = [math]::Round(100.0 * $rows / 9, 1)
  parent_pid = $pidValue
  running = [bool]$process
}
