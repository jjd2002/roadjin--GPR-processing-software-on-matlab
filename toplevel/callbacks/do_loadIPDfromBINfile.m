function [IPD, OPD] = do_loadIPDfromBINfile
%
% Callback function to load the IPD structure from an MGP or MAT file
%
% Copyright (C) 2005,2008 Andreas Tzanis. All rights reserved.
%

global ENVAR
OPD = discardprocdata; 
[BINfname, BINpname]= uigetfile(...
    {'*.mat; *.mgp','MATGPR Files (*.mat, *.mgp)';
     '*.mat','MAT-files (*.mat)'; ...
     '*.mgp', 'MGP-files (*.mgp)'}, ...
    ['Give MATGPR binary data file name'], ...
    ENVAR.currentworkdir);
if BINfname == 0,      % operation canceled
    IPD = initdatastr;
    datafig = findobj('tag','datafigure');
    if ~isempty(datafig),
        delete(datafig);
    end
    return;
end;
if ~isempty(findstr(BINfname,'.mat')),
    IPD = importdata([BINpname BINfname]); 
end
if ~isempty(findstr(BINfname, '.mgp')),
    IPD = readMGP([BINpname BINfname]);
end
% Trap a possible error
if ~isstruct(IPD) || ~ (isfield(IPD,'origin') && isfield(IPD,'sigpos') ...
        && isfield(IPD,'markertr')), 
    erh = errordlg(['The content of this file is incompatible with matGPR!'],...
        'MATGPR : ERROR'); 
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
ENVAR.currentworkdir = BINpname;

return





