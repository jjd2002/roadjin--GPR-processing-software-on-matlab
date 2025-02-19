function OPD = do_splitstepmig(IPD, v2d, zv, dzmig)
%
% Callback function to drive the routine "splitstepmig.m"
% Split-step Phase-shift Depth migration for structures with moderate
% lateral  velocity variations
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
if isempty(v2d); 
    erh = warndlg('Please provide a 2-D velocity model!', ...
        'MATGPR : WARNING!'); 
    uiwait(erh); 
    return; 
end; 
% Proceed
OPD = IPD;
[OPD.d,OPD.z,OPD.ns] = splitstepmig(IPD.d,IPD.dt,IPD.dx,v2d,zv,dzmig);
if isempty(OPD.d),
    disp('SPLITSTEPMIG > Operation aborted - No O/P data returned!');
    return
end
OPD.dz = dzmig; 
OPD.zlab = 'Depth (m)'; 
viewdata(OPD.x,OPD.z,OPD.d,'outdata',OPD.xlab,OPD.zlab); 
title('Split-Step Depth Migrated Section')
iss = size(OPD.history,1); 
OPD.history(iss+1,1) = cellstr('Performed Split-step Depth Migration'); 
return
