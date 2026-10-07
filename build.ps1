param(
    [string]$Compiler = "g++",
    [ValidateSet("Release", "Debug")]
    [string]$Configuration = "Release"
)

$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$binDir = Join-Path $projectRoot "bin"
New-Item -ItemType Directory -Path $binDir -Force | Out-Null

$flags = @("-std=c++17", "-Wall", "-Wextra", "-Wpedantic", "-I", (Join-Path $projectRoot "include"), "-I", (Join-Path $projectRoot "apps"))
if ($Configuration -eq "Release") {
    $flags += @("-O3", "-DNDEBUG")
} else {
    $flags += @("-O0", "-g")
}

$common = @(
    (Join-Path $projectRoot "src\graph.cpp"),
    (Join-Path $projectRoot "src\evaluator.cpp"),
    (Join-Path $projectRoot "src\algorithms.cpp")
)

function Build-Target([string]$Source, [string]$Output) {
    & $Compiler @flags @common $Source "-o" $Output
    if ($LASTEXITCODE -ne 0) {
        throw "Compilation failed for $Source"
    }
}

Build-Target (Join-Path $projectRoot "apps\greedy_main.cpp") (Join-Path $binDir "kmdrdp_greedy.exe")
Build-Target (Join-Path $projectRoot "apps\ga_main.cpp") (Join-Path $binDir "kmdrdp_ga.exe")
Build-Target (Join-Path $projectRoot "apps\vns_main.cpp") (Join-Path $binDir "kmdrdp_vns.exe")
Build-Target (Join-Path $projectRoot "tests\test_core.cpp") (Join-Path $binDir "test_core.exe")

Write-Host "Built binaries in $binDir"

