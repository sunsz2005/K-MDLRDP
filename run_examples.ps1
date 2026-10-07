param([string]$Python = "python")

$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$sample = Join-Path $projectRoot "Instances\sample-2layer"
& (Join-Path $projectRoot "build.ps1")

& (Join-Path $projectRoot "bin\kmdrdp_greedy.exe") $sample 2 13579 10
& (Join-Path $projectRoot "bin\kmdrdp_ga.exe") $sample 2 20 20 2 10 TOURNAMENT 1 13579 1 5 --init greedy
& (Join-Path $projectRoot "bin\kmdrdp_vns.exe") vns $sample 2 13579 1 4 0.4 1 100 50 5 --init greedy
& $Python (Join-Path $projectRoot "python\exact_ilp.py") $sample 2 --model direct --time-limit 30
& $Python (Join-Path $projectRoot "python\exact_ilp.py") $sample 2 --model cumulative --time-limit 30
& $Python (Join-Path $projectRoot "python\ga_ilp_hybrid.py") $sample 2 --time-limit 3 --ga-time 1 --subproblem-time 1 --max-subproblems 10

