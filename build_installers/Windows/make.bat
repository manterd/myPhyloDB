### This file will not install myPhyloDB.
### Please download the installer at: http://www.ars.usda.gov/services/software/download.htm?softwareid=472

cd C:\Users\daniel.manter\Documents\GitHub\myPhyloDB
pyinstaller -D serve-win.spec

cd C:\Users\daniel.manter\Documents\GitHub\myPhyloDB\Installers\Windows
"C:\Program Files (x86)\Inno Setup 5\ISCC.exe" "phyloDB.iss"