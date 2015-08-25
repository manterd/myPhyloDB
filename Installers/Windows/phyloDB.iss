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
Source: "..\..\dist\myPhyloDB\*.pyd"; DestDir: "{app}"; Components: Main
Source: "..\..\dist\myPhyloDB\*.dll"; DestDir: "{app}"; Components: Main
Source: "..\..\dist\myPhyloDB\*.manifest"; DestDir: "{app}"; Components: Main
Source: "..\..\dist\myPhyloDB\*.exe"; DestDir: "{app}"; Components: Main
Source: "..\..\dist\myPhyloDB\dbMicrobe"; DestDir: "{app}"; Flags: uninsneveruninstall; Components: Database
Source: "..\..\dist\myPhyloDB\_MEI\*"; DestDir: "{app}\_MEI"; Flags: recursesubdirs; Components: Main
Source: "..\..\dist\myPhyloDB\django\*"; DestDir: "{app}\django"; Flags: recursesubdirs; Components: Main
Source: "..\..\dist\myPhyloDB\eggs\*"; DestDir: "{app}\eggs"; Flags: recursesubdirs; Components: Main
Source: "..\..\dist\myPhyloDB\include\*"; DestDir: "{app}\include"; Flags: recursesubdirs; Components: Main
Source: "..\..\dist\myPhyloDB\instructions\*"; DestDir: "{app}\instructions"; Flags: recursesubdirs; Components: Manual
Source: "..\..\dist\myPhyloDB\IPython\*"; DestDir: "{app}\IPython"; Flags: recursesubdirs; Components: Main
Source: "..\..\dist\myPhyloDB\media\*"; DestDir: "{app}\media"; Flags: recursesubdirs; Components: Main
Source: "..\..\dist\myPhyloDB\mothur\mothur-win\*"; DestDir: "{app}\mothur\mothur-win"; Flags: recursesubdirs; Components: Main
Source: "..\..\dist\myPhyloDB\mothur\reference\align\*"; DestDir: "{app}\mothur\reference\align"; Flags: recursesubdirs; Components: Main
Source: "..\..\dist\myPhyloDB\mothur\reference\taxonomy\*"; DestDir: "{app}\mothur\reference\taxonomy"; Flags: recursesubdirs; Components: Main
Source: "..\..\dist\myPhyloDB\mothur\reference\template\*"; DestDir: "{app}\mothur\reference\template"; Flags: recursesubdirs; Components: Main
Source: "..\..\dist\myPhyloDB\mpl-data\*"; DestDir: "{app}\mpl-data"; Flags: recursesubdirs; Components: Main
Source: "..\..\dist\myPhyloDB\pytz\*"; DestDir: "{app}\pytz"; Flags: recursesubdirs; Components: Main
Source: "..\..\dist\myPhyloDB\sample_files\*"; DestDir: "{app}\sample_files"; Flags: recursesubdirs; Components: Samples
Source: "..\..\dist\myPhyloDB\sphinx\*"; DestDir: "{app}\sphinx"; Flags: recursesubdirs; Components: Main
Source: "..\..\dist\myPhyloDB\templates\*"; DestDir: "{app}\templates"; Flags: recursesubdirs; Components: Main
Source: "..\..\dist\myPhyloDB\zmq\*"; DestDir: "{app}\zmq"; Flags: recursesubdirs; Components: Main
Source: "..\..\dist\myPhyloDB\R\R-Portable\*"; DestDir: "{app}\R\R-Portable"; Flags: recursesubdirs; Components: Main
Source: "..\..\dist\myPhyloDB\uploads\*"; DestDir: "{app}\uploads"; Flags: recursesubdirs uninsneveruninstall; Components: Uploads

[Icons]
Name: "{group}\myPhyloDB"; Filename: "{app}\myPhyloDB.exe"; IconFilename: "{app}\media\images\database_2_48.ico"
Name: "{group}\Manual"; Filename: "{app}\instructions\Manual.pdf"
Name: "{group}\Uninstall"; Filename: "{uninstallexe}"

[Components]
Name: Main; Description: Core program files
Name: Manual; Description: Instruction manual
Name: Database; Description: Default database
Name: Uploads; Description: Default uploaded files
Name: Samples; Description: Sample files

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



