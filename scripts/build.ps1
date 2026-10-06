$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $PSScriptRoot
$bundledPython = Join-Path $env:USERPROFILE ".cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe"

if (Get-Command python -ErrorAction SilentlyContinue) {
    $pythonExecutable = (Get-Command python).Source
} elseif (Test-Path -LiteralPath $bundledPython) {
    $pythonExecutable = $bundledPython
} else {
    throw "Python 3.11 or newer could not be found."
}

Push-Location $projectRoot
try {
    & $pythonExecutable -m pip wheel . --no-deps --no-build-isolation --wheel-dir dist
    if ($LASTEXITCODE -ne 0) { throw "Wheel build failed." }
    Copy-Item -LiteralPath "examples\two_case.rplab.json" -Destination "dist\two_case.rplab.json" -Force
} finally {
    Pop-Location
}
