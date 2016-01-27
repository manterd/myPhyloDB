rem This file will not install myPhyloDB.
rem To install, please download the installer at: http://www.ars.usda.gov/services/software/download.htm?softwareid=472

cd C:\Users\daniel.manter.NPA-VOICE\PycharmProjects\myPhyloDB
pyinstaller -D serve-win.spec

cd C:\Users\daniel.manter.NPA-VOICE\PycharmProjects\myPhyloDB\build_installers\Windows
"C:\Program Files (x86)\Inno Setup 5\ISCC.exe" "phyloDB.iss"