$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $PSScriptRoot
$instanceRoot = Join-Path (Split-Path -Parent $projectRoot) "k-MLRDP-main\Instances"
$output = Join-Path $projectRoot "Results\formal\n1000_formal_ga_g_hybrid_a.csv"
$seeds = "24680,86420,271828,161803,987654"

Set-Location $projectRoot
& "python" .\python\run_experiments.py $instanceRoot `
    --include-regex "[\\/]n1000(?:-|$)" `
    --methods "ga-g-hybrid" `
    --k protocol `
    --seeds $seeds `
    --time-limit 3600 `
    --output $output
