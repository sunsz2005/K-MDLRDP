$ErrorActionPreference = "Stop"

$root = Split-Path -Parent $PSScriptRoot
$reference = Join-Path (Split-Path -Parent $root) "k-MLRDP-main"
$python = "python"
$executable = Join-Path $root "bin\kmdrdp_vns.exe"
$instance = Join-Path $root "Instances\PPI\PPI-5Layers"
$runsRoot = Join-Path $root "Results\formal\ordinary_roman_vns_certified"
$seeds = 20261006..20261015
$maxParallel = 2

New-Item -ItemType Directory -Force -Path $runsRoot | Out-Null
$pending = [System.Collections.Generic.Queue[int]]::new()
foreach ($seed in $seeds) {
    $runDir = Join-Path $runsRoot "seed_$seed"
    $statusFile = Join-Path $runDir "status.json"
    if (Test-Path -LiteralPath $statusFile) {
        $status = Get-Content -LiteralPath $statusFile -Raw | ConvertFrom-Json
        if ($status.exit_code -eq 0) { continue }
    }
    $pending.Enqueue($seed)
}

$active = @()
while ($pending.Count -gt 0 -or $active.Count -gt 0) {
    while ($pending.Count -gt 0 -and $active.Count -lt $maxParallel) {
        $seed = $pending.Dequeue()
        $runDir = Join-Path $runsRoot "seed_$seed"
        New-Item -ItemType Directory -Force -Path $runDir | Out-Null
        if ($seed -eq 20261006) {
            $greedyCsv = Join-Path $reference "Exec\Greedy\5Layers_seed20261006_k5.init"
        } else {
            $greedyCsv = Join-Path $reference "Results\formal\ordinary_roman_vns_runs\seed_$seed\Greedy\5Layers_seed${seed}_k5.init"
        }
        $initial = Join-Path $runDir "initial.labels"
        & $python (Join-Path $root "python\roman_initial_from_csv.py") $greedyCsv $initial
        if ($LASTEXITCODE -ne 0) { throw "initial-label conversion failed for seed $seed" }

        $stdout = Join-Path $runDir "stdout.log"
        $stderr = Join-Path $runDir "stderr.log"
        $started = Get-Date
        $arguments = @(
            "vns", $instance, "5", "$seed", "2", "12", "0.4", "1",
            "10000", "5000", "3600", "--mode", "roman",
            "--initial-solution", $initial
        )
        $process = Start-Process -FilePath $executable -ArgumentList $arguments `
            -WorkingDirectory $root -WindowStyle Hidden `
            -RedirectStandardOutput $stdout -RedirectStandardError $stderr -PassThru
        $active += [pscustomobject]@{
            seed = $seed; process = $process; started = $started; run_dir = $runDir
        }
    }

    Start-Sleep -Seconds 10
    $stillActive = @()
    foreach ($item in $active) {
        if (-not $item.process.HasExited) {
            $stillActive += $item
            continue
        }
        $finished = Get-Date
        [ordered]@{
            seed = $item.seed
            started = $item.started.ToString("o")
            finished = $finished.ToString("o")
            elapsed_seconds = ($finished - $item.started).TotalSeconds
            exit_code = $item.process.ExitCode
            stdout = Join-Path $item.run_dir "stdout.log"
            stderr = Join-Path $item.run_dir "stderr.log"
        } | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $item.run_dir "status.json") -Encoding UTF8
    }
    $active = $stillActive
}

& $python (Join-Path $root "python\collect_certified_roman_vns.py")
if ($LASTEXITCODE -ne 0) { throw "certificate collection or validation failed" }
$combined = Join-Path $root "Results\formal\ppi_5layer_k5_roman_vns_certified_10seeds.csv"
& $python (Join-Path $root "python\prepare_roman_ppi_consensus.py") $combined `
    --source "certified ordinary-Roman VNS-G-FI, ten independent seeds"
if ($LASTEXITCODE -ne 0) { throw "ordinary-Roman consensus analysis failed" }
& $python (Join-Path $root "python\analyze_roman_double_comparison.py")
if ($LASTEXITCODE -ne 0) { throw "ordinary/double-Roman comparison failed" }
