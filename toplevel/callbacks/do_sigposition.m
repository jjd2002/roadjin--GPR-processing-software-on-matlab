function OPD = do_sigposition(IPD);
%
% Callback function to drive the routine "sigposition.m"
% Adjust signal for time-zero or trim front of time window.
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%

OPD = discardprocdata; 
% Trap common errors
if isempty(IPD.d), 
   erh = errordlg('没有数据可以处理!', 'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end
% Proceed ...
OPD = IPD; 
[OPD.d, OPD.ns, OPD.tt2w] = sigposition('begin',IPD.dt,IPD.dx,IPD.d); 
if isempty(OPD.d),
    disp('SIGPOSITION > Operation aborted - No O/P data returned!');
    return
end
% Operation successful - update signal position info
OPD.sigpos = 0;
% If DZT data, reset header gain function
if ~isempty(IPD.DZThdgain),
    OPD.DZThdgain = IPD.DZThdgain(IPD.ns - OPD.ns + 1 :IPD.ns); 
end
% Display O/P data
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab);
text = ['Moved Time-zero by ' num2str((IPD.ns-OPD.ns)*IPD.dt,4) ' ns'];
title(text);
iss = size(OPD.history,1); 
OPD.history(iss+1,1) = cellstr(text); 
return
