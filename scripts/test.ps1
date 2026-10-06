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

$env:PYTHONPATH = Join-Path $projectRoot "src"
& $pythonExecutable -m unittest discover -s (Join-Path $projectRoot "tests") -v

