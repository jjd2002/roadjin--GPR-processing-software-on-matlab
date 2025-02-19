function IPD = do_saveIPDtoMATfile(IPD)
%
% Callback function to save the IPD structure to a MAT file
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
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
MATfname = [IPD.fname(1:getts) num2str(IPD.TimesSaved) '.mat'];
[MATfname, MATpname]= uiputfile('*.mat; *.MAT', ...
    'Save Current Input Data as ...',[ENVAR.currentworkdir MATfname]); 
if MATfname == 0;          % Cancel has been pressed
    IPD.TimesSaved = IPD.TimesSaved - 1; 
    return; 
end; 
IPD.fname = MATfname; 
IPD.pname = MATpname;
save([MATpname MATfname],'IPD',ENVAR.matfileformat); 
return
