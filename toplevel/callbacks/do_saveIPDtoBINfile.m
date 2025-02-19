function do_saveIPDtoBINfile(IPD)
%
% Callback function to save the IPD structure to an MGP or MAT file
%
% Copyright (C) 2005,2008 Andreas Tzanis. All rights reserved.
%

global ENVAR
% Trap common errors
if isempty(IPD.d), 
   erh = errordlg('No data in the workspace to save!', 'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end
% Proceed ...
IPD.TimesSaved = IPD.TimesSaved + 1; 
getts  = findstr(IPD.fname,'.');
BINfname = [IPD.fname(1:getts) num2str(IPD.TimesSaved) '.' ENVAR.binary_file_type];
if strcmp(ENVAR.binary_file_type,'mgp'),
    Ofiles = {'*.mgp;', 'MGP-files (*.mgp)';
     '*.mat','MAT-files (*.mat)'; ...
     '*.mgp', 'MGP-files (*.mgp)'};
elseif strcmp(ENVAR.binary_file_type,'mat'),
    Ofiles = {'*.mat;', 'MAT-files (*.mat)';
     '*.mgp', 'MGP-files (*.mgp)'; ...
     '*.mat','MAT-files (*.mat)'};
end
[BINfname, BINpname]= uiputfile( Ofiles, ...
    'Save Current Input Data as ...',[ENVAR.currentworkdir BINfname]); 
if BINfname == 0;          % Cancel has been pressed
    IPD.TimesSaved = IPD.TimesSaved - 1; 
    return; 
end; 
IPD.fname = BINfname; 
IPD.pname = BINpname;
if strcmp(ENVAR.binary_file_type, 'mat'),
    save([BINpname BINfname],'IPD',ENVAR.matfileformat); 
elseif strcmp(ENVAR.binary_file_type, 'mgp'),
    save2MGP(IPD,[BINpname BINfname]);
end
return
