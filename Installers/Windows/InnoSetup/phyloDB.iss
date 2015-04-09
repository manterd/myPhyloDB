[Setup]
AppName=myPhyloDB
AppVersion=1.0
DefaultDirName={code:DefDirRoot}\myPhyloDB
DefaultGroupName=myPhyloDB
Compression=lzma2
SolidCompression=yes
ArchitecturesAllowed=x64
ArchitecturesInstallIn64BitMode=x64
PrivilegesRequired=lowest
OutputBaseFilename=myPhyloDB_1.0_Win_x64_install

[Files]
Source: "..\..\..\dist\myPhyloDB\*"; DestDir: "{app}";
Source: "..\..\..\dist\myPhyloDB\_MEI\*"; DestDir: "{app}\_MEI"; Flags: recursesubdirs
Source: "..\..\..\dist\myPhyloDB\django\*"; DestDir: "{app}\django"; Flags: recursesubdirs
Source: "..\..\..\dist\myPhyloDB\eggs\*"; DestDir: "{app}\eggs"; Flags: recursesubdirs
Source: "..\..\..\dist\myPhyloDB\include\*"; DestDir: "{app}\include"; Flags: recursesubdirs
Source: "..\..\..\dist\myPhyloDB\instructions\*"; DestDir: "{app}\instructions"; Flags: recursesubdirs
Source: "..\..\..\dist\myPhyloDB\IPython\*"; DestDir: "{app}\IPython"; Flags: recursesubdirs
Source: "..\..\..\dist\myPhyloDB\media\*"; DestDir: "{app}\media"; Flags: recursesubdirs
Source: "..\..\..\dist\myPhyloDB\mpl-data\*"; DestDir: "{app}\mpl-data"; Flags: recursesubdirs
Source: "..\..\..\dist\myPhyloDB\pytz\*"; DestDir: "{app}\pytz"; Flags: recursesubdirs
Source: "..\..\..\dist\myPhyloDB\sample_files\*"; DestDir: "{app}\sample_files"; Flags: recursesubdirs
Source: "..\..\..\dist\myPhyloDB\sphinx\*"; DestDir: "{app}\sphinx"; Flags: recursesubdirs
Source: "..\..\..\dist\myPhyloDB\templates\*"; DestDir: "{app}\templates"; Flags: recursesubdirs
Source: "..\..\..\dist\myPhyloDB\zmq\*"; DestDir: "{app}\zmq"; Flags: recursesubdirs

[Icons]
Name: "{group}\myPhyloDB"; Filename: "{app}\myPhyloDB.exe"; IconFilename: "{app}\media\images\database_2_48.ico"
Name: "{group}\Manual"; Filename: "{app}\instructions\Manual.pdf"
Name: "{group}\Uninstall"; Filename: "{uninstallexe}"

[Code]
function IsRegularUser(): Boolean;
begin
Result := not (IsAdminLoggedOn or IsPowerUserLoggedOn);
end;

function DefDirRoot(Param: String): String;
begin
if IsRegularUser then
Result := ExpandConstant('{localappdata}')
else
Result := ExpandConstant('{pf}')
end;



