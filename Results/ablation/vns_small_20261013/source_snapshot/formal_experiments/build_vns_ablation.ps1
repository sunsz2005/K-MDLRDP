$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $PSScriptRoot
Set-Location $projectRoot
$common = @("src/graph.cpp", "src/evaluator.cpp", "src/algorithms.cpp")
& g++ -std=c++17 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Iinclude -Iapps @common apps/vns_ablation_main.cpp -o bin/vns_ablation.exe
if ($LASTEXITCODE -ne 0) { throw "ablation compilation failed" }
& g++ -std=c++17 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Iinclude src/graph.cpp src/evaluator.cpp apps/evaluation_benchmark_main.cpp -o bin/evaluation_benchmark.exe
if ($LASTEXITCODE -ne 0) { throw "benchmark compilation failed" }
& g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -Iinclude @common tests/test_ablation.cpp -o bin/test_ablation.exe
if ($LASTEXITCODE -ne 0) { throw "ablation tests compilation failed" }
& g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -Iinclude @common tests/test_core.cpp -o bin/test_core_ablation.exe
if ($LASTEXITCODE -ne 0) { throw "core tests compilation failed" }
& ./bin/test_ablation.exe
if ($LASTEXITCODE -ne 0) { throw "ablation tests failed" }
& ./bin/test_core_ablation.exe
if ($LASTEXITCODE -ne 0) { throw "core tests failed" }
