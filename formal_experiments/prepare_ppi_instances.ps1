$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $PSScriptRoot
$sourceRoot = Join-Path (Split-Path -Parent $projectRoot) "k-MLRDP-main\Instances\PPI"
$targetRoot = Join-Path $projectRoot "Instances\PPI"
$reportPath = Join-Path $projectRoot "data_reconstruction\output\reconstruction_report.json"

$report = Get-Content -LiteralPath $reportPath -Raw | ConvertFrom-Json
foreach ($count in 3, 4, 5) {
    $status = $report.mapping_status."${count}Layers"
    if ($status -ne "complete and exact") {
        throw "PPI ${count}-layer mapping is not complete and exact: $status"
    }
}

$layers = @{
    3 = @("STRING", "iRefIndex", "wiphi")
    4 = @("STRING", "iRefIndex", "wiphi", "BioGRID")
    5 = @("STRING", "iRefIndex", "wiphi", "BioGRID", "MINT")
}

New-Item -ItemType Directory -Force -Path $targetRoot | Out-Null
$manifestInstances = @()
foreach ($count in 3, 4, 5) {
    $sourceDir = Join-Path $sourceRoot "${count}Layers"
    $targetDir = Join-Path $targetRoot "PPI-${count}Layers"
    New-Item -ItemType Directory -Force -Path $targetDir | Out-Null
    for ($index = 0; $index -lt $layers[$count].Count; $index++) {
        $source = Join-Path $sourceDir ($layers[$count][$index] + "_mapped_names.txt")
        $target = Join-Path $targetDir ("layer_{0}.txt" -f ($index + 1))
        if (-not (Test-Path -LiteralPath $source)) {
            throw "Missing source layer: $source"
        }
        if (-not (Test-Path -LiteralPath $target)) {
            New-Item -ItemType HardLink -Path $target -Target $source | Out-Null
        }
        if ((Get-Item -LiteralPath $source).Length -ne (Get-Item -LiteralPath $target).Length) {
            throw "Layer length mismatch: $target"
        }
    }
    $manifestLayers = @()
    for ($index = 0; $index -lt $layers[$count].Count; $index++) {
        $target = Join-Path $targetDir ("layer_{0}.txt" -f ($index + 1))
        $header = (Get-Content -LiteralPath $target -TotalCount 1) -split "\s+"
        $manifestLayers += [ordered]@{
            index = $index + 1
            source = $layers[$count][$index]
            vertices = [int]$header[0]
            edges = [int]$header[1]
            sha256 = (Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash
        }
    }
    $manifestInstances += [ordered]@{
        instance = "PPI-${count}Layers"
        mapping_status = $report.mapping_status."${count}Layers"
        layers = $manifestLayers
    }
}

$manifest = [ordered]@{
    generated_on = (Get-Date -Format "yyyy-MM-dd")
    reconstruction_report = "data_reconstruction/output/reconstruction_report.json"
    instances = $manifestInstances
}
$manifest | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (
    Join-Path $targetRoot "manifest.json"
) -Encoding UTF8

Get-ChildItem -LiteralPath $targetRoot -Directory | ForEach-Object {
    [pscustomobject]@{
        Instance = $_.Name
        Layers = @(Get-ChildItem -LiteralPath $_.FullName -Filter "layer_*.txt").Count
        Bytes = (Get-ChildItem -LiteralPath $_.FullName -Filter "layer_*.txt" |
            Measure-Object -Property Length -Sum).Sum
    }
}
