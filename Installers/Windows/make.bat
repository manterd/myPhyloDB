cd C:\Users\daniel.manter\Documents\GitHub\myPhyloDB
pyinstaller -D serve-win.spec

cd C:\Users\daniel.manter\Documents\GitHub\myPhyloDB\Installers\Windows
"C:\Program Files (x86)\Inno Setup 5\ISCC.exe" "phyloDB.iss"