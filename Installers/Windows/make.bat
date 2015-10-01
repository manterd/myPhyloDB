cd C:\Users\daniel.manter.NPA-VOICE\PycharmProjects\myPhyloDB
pyinstaller -D serve-win.spec

cd C:\Users\daniel.manter.NPA-VOICE\PycharmProjects\myPhyloDB\Installers\Windows
"C:\Program Files (x86)\Inno Setup 5\ISCC.exe" "phyloDB.iss"