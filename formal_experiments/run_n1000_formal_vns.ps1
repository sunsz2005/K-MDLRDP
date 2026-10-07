$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $PSScriptRoot
$instanceRoot = Join-Path (Split-Path -Parent $projectRoot) "k-MLRDP-main\Instances"
$output = Join-Path $projectRoot "Results\formal\n1000_formal_vns_g_fi_extra.csv"
$seeds = "24680,97531,86420,112358,271828,314159,161803,424242,987654"

Set-Location $projectRoot
& "python" .\python\run_experiments.py $instanceRoot `
    --include-regex "[\\/]n1000(?:-|$)" `
    --methods "vns-g-fi" `
    --k protocol `
    --seeds $seeds `
    --time-limit 3600 `
    --output $output
