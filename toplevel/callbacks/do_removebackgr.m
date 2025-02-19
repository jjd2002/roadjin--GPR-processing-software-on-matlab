function OPD = do_removebackgr(IPD);
%
% Callback function to drive global background trace removal 
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%

if ~isempty(IPD.d),
    OPD = discardprocdata; 
    OPD = IPD;
    OPD.d = rmbackgr(IPD.d); 
    if ~isempty(OPD.d),
        viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab); 
        iss = size(OPD.history,1); 
        OPD.history(iss+1,1) = cellstr('Removed Global Background Trace'); 
    end;
else
    erh = errordlg('No data to process!','MATGPR: ERROR');
    uiwait(erh)
    return
end
return
