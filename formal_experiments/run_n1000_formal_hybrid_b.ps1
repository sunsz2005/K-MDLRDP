$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $PSScriptRoot
$instanceRoot = Join-Path (Split-Path -Parent $projectRoot) "k-MLRDP-main\Instances"
$output = Join-Path $projectRoot "Results\formal\n1000_formal_ga_g_hybrid_b.csv"
$seeds = "97531,112358,314159,424242"

Set-Location $projectRoot
& "python" .\python\run_experiments.py $instanceRoot `
    --include-regex "[\\/]n1000(?:-|$)" `
    --methods "ga-g-hybrid" `
    --k protocol `
    --seeds $seeds `
    --time-limit 3600 `
    --output $output
