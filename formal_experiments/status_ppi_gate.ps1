$projectRoot = Split-Path -Parent $PSScriptRoot
$resultRoot = Join-Path $projectRoot "Results\formal"
$jobs = @(
    @{ Name = "Greedy"; File = "ppi_gate_greedy.csv"; Expected = 12; Pid = "ppi_gate_greedy.pid" },
    @{ Name = "VNS-G-FI"; File = "ppi_gate_vns_g_fi.csv"; Expected = 12; Pid = "ppi_gate_vns.pid" },
    @{ Name = "GA-G-Hybrid"; File = "ppi_gate_ga_g_hybrid.csv"; Expected = 12; Pid = "ppi_gate_hybrid.pid" }
)

foreach ($job in $jobs) {
    $resultPath = Join-Path $resultRoot $job.File
    $pidPath = Join-Path $resultRoot $job.Pid
    $rows = if (Test-Path -LiteralPath $resultPath) { @(Import-Csv -LiteralPath $resultPath).Count } else { 0 }
    $processId = if (Test-Path -LiteralPath $pidPath) { [int](Get-Content -LiteralPath $pidPath -Raw) } else { 0 }
    $running = $processId -gt 0 -and $null -ne (Get-Process -Id $processId -ErrorAction SilentlyContinue)
    [pscustomobject]@{
        Job = $job.Name
        Rows = $rows
        Expected = $job.Expected
        Percent = [math]::Round(100 * $rows / $job.Expected, 1)
        ProcessId = $processId
        Running = $running
    }
}
