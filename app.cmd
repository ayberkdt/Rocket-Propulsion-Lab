@echo off
rem Rocket Propulsion Lab - cift tiklayinca yerel hesap sunucusunu baslatir
rem ve calisma tezgahini varsayilan tarayicide acar. Pencereyi kapatmak
rem sunucuyu durdurur. Zaten calisiyorsa yalnizca tarayiciyi acar.
setlocal
chcp 65001 >nul
cd /d "%~dp0"
set "PYTHONPATH=%~dp0src"
title Rocket Propulsion Lab

set "PYTHON="
where py >nul 2>nul && py -3 -c "import sys; sys.exit(sys.version_info < (3, 11))" >nul 2>nul && set "PYTHON=py -3"
if not defined PYTHON (
  where python >nul 2>nul && python -c "import sys; sys.exit(sys.version_info < (3, 11))" >nul 2>nul && set "PYTHON=python"
)
if not defined PYTHON (
  echo.
  echo Python 3.11 veya daha yenisi bulunamadi.
  echo https://www.python.org/downloads/ adresinden kurup bu dosyayi yeniden calistirin.
  echo.
  pause
  exit /b 1
)

echo Rocket Propulsion Lab baslatiliyor...
%PYTHON% -m rocket_propulsion.api.server --open %*

if errorlevel 1 (
  echo.
  echo Sunucu beklenmedik sekilde durdu. Yukaridaki mesaji kontrol edin.
  pause
)
endlocal
