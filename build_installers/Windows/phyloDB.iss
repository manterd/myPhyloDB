[Setup]
AppName=myPhyloDB
AppVersion=1.2.0
DefaultDirName={code:DefDirRoot}\myPhyloDB
DefaultGroupName=myPhyloDB
Compression=lzma2
SolidCompression=yes
ArchitecturesAllowed=x64
ArchitecturesInstallIn64BitMode=x64
PrivilegesRequired=lowest
OutputBaseFilename=myPhyloDB_v.1.2.0_Win_x64_install


[Files]
Source: "..\..\dist\myPhyloDB\db.Microbe"; DestDir: "{app}"; Flags: uninsneveruninstall; Components: Database
Source: "..\..\dist\myPhyloDB\db.PICRUSt"; DestDir: "{app}"; Flags: uninsneveruninstall; Components: Database

Source: "..\..\dist\myPhyloDB\*.pyd"; DestDir: "{app}"; Components: Main
Source: "..\..\dist\myPhyloDB\*.dll"; DestDir: "{app}"; Components: Main
Source: "..\..\dist\myPhyloDB\*.manifest"; DestDir: "{app}"; Components: Main
Source: "..\..\dist\myPhyloDB\*.exe"; DestDir: "{app}"; Components: Main

Source: "..\..\dist\myPhyloDB\Include\*"; DestDir: "{app}\Include"; Flags: recursesubdirs; Components: Main
Source: "..\..\dist\myPhyloDB\R\R-portable\*"; DestDir: "{app}\R\R-Portable"; Flags: recursesubdirs; Components: Main
Source: "..\..\dist\myPhyloDB\config\*"; DestDir: "{app}\config"; Flags: recursesubdirs; Components: Main
Source: "..\..\dist\myPhyloDB\django\*"; DestDir: "{app}\django"; Flags: recursesubdirs; Components: Main
Source: "..\..\dist\myPhyloDB\instructions\current\*"; DestDir: "{app}\instructions\current"; Flags: recursesubdirs; Components: Manual
Source: "..\..\dist\myPhyloDB\media\*"; DestDir: "{app}\media"; Flags: recursesubdirs; Components: Main
Source: "..\..\dist\myPhyloDB\mothur\mothur-win\*"; DestDir: "{app}\mothur\mothur-win"; Flags: recursesubdirs; Components: Main
Source: "..\..\dist\myPhyloDB\mothur\reference\zip\*"; DestDir: "{app}\mothur\reference\zip"; Flags: recursesubdirs; Components: Main
Source: "..\..\dist\myPhyloDB\pytz\*"; DestDir: "{app}\pytz"; Flags: recursesubdirs; Components: Main
Source: "..\..\dist\myPhyloDB\templates\*"; DestDir: "{app}\templates"; Flags: recursesubdirs; Components: Main
Source: "..\..\dist\myPhyloDB\uploads\*"; DestDir: "{app}\uploads"; Flags: recursesubdirs uninsneveruninstall; Components: Main


[Icons]
Name: "{group}\myPhyloDB"; Filename: "{app}\myPhyloDB.exe"; IconFilename: "{app}\media\images\myPhyloDB_Logo.ico"
Name: "{group}\Manual.pdf"; Filename: "{app}\instructions\current\Manual.pdf"
Name: "{group}\Uninstall"; Filename: "{uninstallexe}"


[Components]
Name: "Main"; Description: "Core program files"; Types: full compact;
Name: "Database"; Description: "Default database"; Types: full;
Name: "Manual"; Description: "Help Manual"; Types: full;


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



