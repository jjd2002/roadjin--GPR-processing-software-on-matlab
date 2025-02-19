function OPD = do_pspimig(IPD, v2d, zv, dzmig)
%
% Callback function to drive the routine "migpspi.m"
% Phase-shift plus Interpolation Depth migration for structures with
% lateral  velocity variations.
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
% Velocity is a 3 element celll array, but pspimig requires only the
% non-dispersive term passed in v2d{1}.
if iscell(v2d),
    [OPD.d, OPD.z, OPD.ns] = pspimig(flipud(IPD.d), IPD.dt, IPD.dx, v2d{1}, zv, dzmig);
else
    [OPD.d, OPD.z, OPD.ns] = pspimig(flipud(IPD.d), IPD.dt, IPD.dx, v2d, zv, dzmig);
end
if isempty(OPD.d),
    disp('MIGPSPI > Operation aborted - No O/P data returned!');
    return
end
OPD.dz   = dzmig; 
OPD.zlab = 'Depth (m)'; 
viewdata(OPD.x,OPD.z,OPD.d,'outdata',OPD.xlab,OPD.zlab); 
title('Phase-shift + Interpolation Depth Migrated Section')
iss = size(OPD.history,1); 
text = 'Performed Phase-shift + Interpolation Depth Migration';
OPD.history(iss+1,1) = cellstr(text); 
return
