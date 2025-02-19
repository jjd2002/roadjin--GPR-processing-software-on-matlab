function OPD = do_trimtimewindow(IPD)
%
% Callback function to drive the routine "trimtimewindow.m"
% Trim the rear end of the recorded time window
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%

OPD = discardprocdata; 
% Trap common errors
if isempty(IPD.d), 
   erh = errordlg('请先导入数据!', 'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end
% Proceed ...
OPD = IPD;
[OPD.d,OPD.tt2w,OPD.ns] = trimtimewindow(IPD.d,IPD.tt2w,IPD.dt,IPD.ns); 
if isempty(OPD.d),
    disp('TRIM TIME WINDOW > Operation aborted - No data returned!');
    return
end
% If DZT data, trim the header gain information as well!
if ~isempty(IPD.DZThdgain), 
    OPD.DZThdgain = OPD.DZThdgain(1:OPD.ns);
end
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab); 
% Update history
iss = size(OPD.history,1); 
text = ['Trimmed Time Window; Discarded last ' ...
    num2str((IPD.ns-OPD.ns)*IPD.dt,5) ' ns (' num2str(IPD.ns-OPD.ns) ...
    ' samples)'];
OPD.history(iss+1,1) = cellstr(text); 
return
