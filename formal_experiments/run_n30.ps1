param([string]$Python = "python")

$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$instanceRoot = Join-Path (Split-Path -Parent $projectRoot) "k-MLRDP-main\Instances"
$runner = Join-Path $projectRoot "python\run_experiments.py"
$output = Join-Path $projectRoot "Results\formal\n30_runs.csv"
$methods = "dr-new1,dr-new2,greedy,vns-r-fi,vns-r-bi,ga-r-fi,ga-r-hybrid,vns-g-fi,vns-g-bi,ga-g-fi,ga-g-hybrid"
$seeds = "13579,24680,97531,86420,112358,271828,314159,161803,424242,987654"

& $Python $runner $instanceRoot `
    --methods $methods `
    --k protocol `
    --seeds $seeds `
    --time-limit 3600 `
    --include-regex 'n30(?:-|$)' `
    --output $output
exit $LASTEXITCODE
