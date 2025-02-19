function [IPD, OPD] = do_loadIPDfromMATfile(IPD)
%
% Callback function to load the IPD structure from a MAT file
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%

global ENVAR
OPD = discardprocdata; 
% Get file 
[MATfname, MATpname]= uigetfile('*.mat; *.MAT', ...
    'Give MAT data file name', ENVAR.currentworkdir);
if MATfname == 0,      % operation canceled
    IPD = initdatastr;
    datafig = findobj('tag','datafigure');
    if ~isempty(datafig),
        delete(datafig);
    end
    return;
end;
IPD = importdata([MATpname MATfname]); 
% Trap a possible error
if isempty(IPD.d) && ~ischar(IPD.fname), 
    erh = errordlg(['The content of this MAT file is incompatible ' ...
        'with matGPR!'], 'MATGPR : ERROR'); 
    uiwait(erh); 
    IPD = initdatastr;
    return; 
end; 
% Project data info on the MATGPR information window
showinfo(IPD); 
set(findobj('tag','undoxx'),'Enable','off');
set(findobj('tag','undoundo'),'Enable','off');

% Erase the UNDO/RESTORE buffer of the previous session, if any
set(findobj('tag','fi0'),'userdata',[]);

% Display the IPD
if isempty(IPD.z),
    viewdata(IPD.x,IPD.tt2w,IPD.d,'indata',IPD.xlab,IPD.zlab); 
else
    viewdata(IPD.x,IPD.z,IPD.d,'indata',IPD.xlab,IPD.zlab); 
end

% Reset the current work directory
ENVAR.currentworkdir = MATpname;

return





