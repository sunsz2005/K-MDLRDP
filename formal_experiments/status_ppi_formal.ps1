$projectRoot = Split-Path -Parent $PSScriptRoot
$resultFile = Join-Path $projectRoot "Results\formal\ppi_formal_ga_g_hybrid_k1_extra.csv"
$pidFile = Join-Path $projectRoot "Results\formal\ppi_formal_hybrid_k1.pid"

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
  expected_new_runs = 27
  percent = [math]::Round(100.0 * $rows / 27, 1)
  parent_pid = $pidValue
  running = [bool]$process
}
