$projectRoot = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$result = Join-Path $projectRoot "Results\formal\n30_runs.csv"
$pidFile = Join-Path $projectRoot "Results\formal\n30.pid"
$expected = 2208

$pidValue = if (Test-Path $pidFile) { (Get-Content $pidFile -Raw).Trim() } else { "" }
$running = $false
if ($pidValue) {
    $running = $null -ne (Get-Process -Id ([int]$pidValue) -ErrorAction SilentlyContinue)
}
$rows = if (Test-Path $result) { @(Import-Csv $result) } else { @() }
$ok = @($rows | Where-Object status -eq "ok").Count
$failed = @($rows | Where-Object status -eq "error").Count

[pscustomobject]@{
    ProcessId = $pidValue
    Running = $running
    CompletedRows = $rows.Count
    SuccessfulRows = $ok
    FailedRows = $failed
    ExpectedRows = $expected
    CompletionPercent = [math]::Round(100.0 * $rows.Count / $expected, 2)
} | Format-List

if ($rows.Count -gt 0) {
    $rows | Group-Object method,status | Sort-Object Name |
        Select-Object Name,Count | Format-Table -AutoSize
}
