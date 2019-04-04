@echo off


:main_fun
    call:global_variables
    if [%1] == [clean] (
        call:clean_kernel
        call:clean_package
        goto:eof
    )

    if [%1] == [kernel] (
        call:compile_kernel
        call:copy_depends
        goto:eof
    ) else (
        call:usage
    )
goto:eof

:global_variables
    set x64=1
    set MINGW_HOME=C:\D\Work\code\software\mingw64
    set MSYS_HOME=C:\D\Work\code\software\mingw64\msys\1.0

    if not exist %MINGW_HOME% echo please install minGW first || pause
    if not exist %MSYS_HOME% echo please install msys 1.0 first || pause

    set PATH=%MINGW_HOME%\bin;%MSYS_HOME%\bin;%DOTNET_HOME%\;%NSIS%\;%TortoiseSVN%\bin;
    set INCLUDE=%MINGW_HOME%\include
    set C_INCLUDE_PATH=%MINGW_HOME%\include
    set CPLUS_INCLUDE_PATH=%MINGW_HOME%\include
    set LIBRARY_PATH=%MINGW_HOME%\lib

    set CUR_DIR=%cd%

    set CONFIGS=%CUR_DIR%\kernel\Trunk\configs\*
    set THIRD_PARTYS=%CUR_DIR%\kernel\Trunk\third_party\*

    
goto:eof

:usage
    echo Usage: make.bat [clean^|kernel]
    pause
goto:eof


:compile_kernel
    cd kernel\Trunk
    if [%x64%] == [0] (
        make BIT64=0 UNITTEST=0 BUILD_DIR=%CUR_DIR:\=/%/build INSTALL_DIR=%CUR_DIR:\=/%/bin CONFIG=%CUR_DIR:\=/%/kernel/Trunk/config
        make install BIT64=0 UNITTEST=0 BUILD_DIR=%CUR_DIR:\=/%/build INSTALL_DIR=%CUR_DIR:\=/%/bin CONFIG=%CUR_DIR:\=/%/kernel/Trunk/config
    ) else (
        make BIT64=1 UNITTEST=0 BUILD_DIR=%CUR_DIR:\=/%/build INSTALL_DIR=%CUR_DIR:\=/%/bin CONFIG=%CUR_DIR:\=/%/kernel/Trunk/config
        make install BIT64=1 UNITTEST=0 BUILD_DIR=%CUR_DIR:\=/%/build INSTALL_DIR=%CUR_DIR:\=/%/bin CONFIG=%CUR_DIR:\=/%/kernel/Trunk/config
    )
    cd %CUR_DIR%
goto:eof


:clean_kernel
    cd kernel\Trunk
    echo begin to clean kernel...
    make clean BIT64=0 UNITTEST=0 BUILD_DIR=%CUR_DIR:\=/%/build INSTALL_DIR=%CUR_DIR:\=/%/bin CONFIG=%CUR_DIR:\=/%/kernel/Trunk/config
    echo finished kernel cleaning.
    cd %CUR_DIR%
goto:eof


:clean_package
    cd installer
    rd /S /Q pSimulate
    cd %CUR_DIR%
goto:eof

:copy_depends
    cp %CONFIGS:\=/% %CUR_DIR%/bin
	cp %THIRD_PARTYS:\=/% %CUR_DIR%/bin
	
    if not exist %CUR_DIR%\installer\pSimulate\bin md %CUR_DIR%\installer\pSimulate\bin
    cp %CUR_DIR:\=/%/bin/* %CUR_DIR:\=/%/installer/pSimulate/bin
	
	if not exist %CUR_DIR%\installer\pSimulate\doc md %CUR_DIR%\installer\pSimulate\doc
	cp %CUR_DIR%\installer\pSimulate\bin\"pSimulate User Guide.pdf" %CUR_DIR:\=/%/installer/pSimulate/doc
	
	cp -rf %CUR_DIR:\=/%/example %CUR_DIR:\=/%/installer/pSimulate
goto:eof

