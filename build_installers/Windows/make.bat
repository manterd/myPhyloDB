rem This file is used to make a new myPhyloDB installer
rem it will not install myPhyloDB

rem If you are looking for a pre-built installer
rem Please download the installer at: http://www.ars.usda.gov/services/software/download.htm?softwareid=472

rem To run this file type the following in your terminal
rem cd C:\Users\daniel.manter.\Documents\GitHub\myPhyloDB
rem make.bat

cd C:\Users\daniel.manter.\Documents\GitHub\myPhyloDB
pyinstaller -D serve-win.spec

copy "C:\Users\daniel.manter\AppData\Local\Continuum\Miniconda2\Library\bin\mkl_avx2.dll" "C:\Users\daniel.manter\Documents\GitHub\myPhyloDB\dist\myPhyloDB\mkl_avx2.dll"
copy "C:\Users\daniel.manter\AppData\Local\Continuum\Miniconda2\Library\bin\mkl_def.dll" "C:\Users\daniel.manter\Documents\GitHub\myPhyloDB\dist\myPhyloDB\mkl_def.dll"

cd C:\Users\daniel.manter.\Documents\GitHub\myPhyloDB\build_installers\Windows
"C:\Program Files (x86)\Inno Setup 5\ISCC.exe" "phyloDB.iss"