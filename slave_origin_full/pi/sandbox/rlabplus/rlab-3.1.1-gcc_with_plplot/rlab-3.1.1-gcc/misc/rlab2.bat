@echo off

REM
REM Batch file for running RLaB in a win32 environment.
REM

REM
REM Set the installation directory.
REM YOU MUST EDIT THE FOLLOWING LINE TO MATCH YOUR INSTALLATION
REM

set INSD=c:\rlab2

REM
REM Setup the runtime environment
REM

set RLAB2_PATH=.;%INSD%\rlib;%INSD%\toolbox;%INSD%\examples
set RLAB2_RC0=%INSD%\rlab.rc
set RLAB2_LIB_DIR=%INSD%\rlib
set RLAB2_HELP_DIR=%INSD%\doc\help
set RLAB2_PAGER="MORE"
set RLAB2_HELP_PAGER="MORE"

REM
REM Setup Plplot environment
REM

set PLPLOT_HOME=%INSD%
set PLPLOT_HOME_ENV=%INSD%
set PLPLOT_BIN_ENV=%INSD%
set PLPLOT_LIB=%INSD%
set PLPLOT_LIB_ENV=%INSD%
set PLPLOT_LIB_DIR=%INSD%
set PLPLOT_DIR=%INSD%
set PLPLOT_FONT_DIR=%INSD%

REM
REM Run rlab...
REM

%INSD%\rlab %1 %2 %3 %4 %5 %6 %7 %8 %9