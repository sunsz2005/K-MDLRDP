$projectRoot = Split-Path -Parent $PSScriptRoot
$resultRoot = Join-Path $projectRoot "Results\formal"
$jobs = @(
    @{ Name = "VNS-G-FI"; File = "n1000_formal_vns_g_fi_extra.csv"; Expected = 216; Pid = "n1000_formal_vns.pid" },
    @{ Name = "GA-G-Hybrid-A"; File = "n1000_formal_ga_g_hybrid_a.csv"; Expected = 120; Pid = "n1000_formal_hybrid_a.pid" },
    @{ Name = "GA-G-Hybrid-B"; File = "n1000_formal_ga_g_hybrid_b.csv"; Expected = 96; Pid = "n1000_formal_hybrid_b.pid" }
)

foreach ($job in $jobs) {
    $resultPath = Join-Path $resultRoot $job.File
    $pidPath = Join-Path $resultRoot $job.Pid
    $rows = if (Test-Path -LiteralPath $resultPath) { @(Import-Csv -LiteralPath $resultPath).Count } else { 0 }
    $processId = if (Test-Path -LiteralPath $pidPath) { [int](Get-Content -LiteralPath $pidPath -Raw) } else { 0 }
    $running = $false
    if ($processId -gt 0) {
        $running = $null -ne (Get-Process -Id $processId -ErrorAction SilentlyContinue)
    }
    [pscustomobject]@{
        Job = $job.Name
        Rows = $rows
        Expected = $job.Expected
        Percent = [math]::Round(100 * $rows / $job.Expected, 1)
        ProcessId = $processId
        Running = $running
    }
}
