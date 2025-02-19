function OPD = do_markerinterp(IPD)
%
% Callback function to drive marker interpolation from unequally to equally
% spaced traces.
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%

global ENVAR
OPD = discardprocdata; 
% This function does not handle bistatic data - Check IPD.origin ...
type1  = findstr(IPD.origin,'Bistatic');
type2  = findstr(IPD.origin,'PULSE EKKO');
thesis = ~isempty(type1) || ~isempty(type2);
if thesis,  
    erh = errordlg('Function presently not available for BISTATIC data',...
        'MATGPR: ERROR');
    uiwait(erh);
    return
end
% Trap common errors
if isempty(IPD.d), 
   erh = errordlg('No data to process!', 'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end
if ~isempty(IPD.dx), 
    erh = msgbox({'The data is already equally-spaced! ' 
        'Use resampling instead! '},'MATGPR: INFO','help');
    uiwait(erh); 
    return;
end
OPD = IPD; 
[OPD.d,OPD.ntr,OPD.dx,OPD.markertr,OPD.x] = ...
    markerinterp(IPD.d, IPD.markertr, ENVAR.currentworkdir, IPD.fname); 
if isempty(OPD.d),
    disp('MARKER INTERPOLATION > Operation aborted - No data returned!');
    return
end
OPD.xlab = 'Scan Axis (meters)'; 
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab); 
iss = size(OPD.history,1); 
text = ['Transformed to Equal Spacing; Traces ' num2str(OPD.ntr) ...
    '; Spacing ' num2str(OPD.dx,4) ' m'];
OPD.history(iss+1,1) = cellstr(text); 

return
