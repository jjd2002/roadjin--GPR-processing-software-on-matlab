function OPD = do_stoltmig(IPD, v1d)
%
% Callback function to drive the routine "stoltmig.m"
% Stolt F-K migration for constant or layered velocity structures
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%

OPD = discardprocdata; 
% Trap common errors ...
if isempty(IPD.d), 
   erh = errordlg('No data to process!', 'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end
if isempty(IPD.dx), 
    erh = errordlg('Operation possible only for equally-spaced data!', ...
        'MATGPR : ERROR');
    uiwait(erh); 
    return;
end; 
if isempty(v1d); 
    erh = warndlg('Please provide an 1-D velocity model!', ...
        'MATGPR : WARNING'); 
    uiwait(erh); 
    return; 
end; 
% Proceed
OPD = IPD;
OPD.d = stoltmig(IPD.d, IPD.dt, IPD.dx, v1d);
if isempty(OPD.d),
    disp('STOLTMIG > Operation aborted - No O/P data returned!');
    return
end
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab); 
iss = size(OPD.history,1); 
OPD.history(iss+1,1) = cellstr('Done F-K (Stolt) Migration'); 
return
