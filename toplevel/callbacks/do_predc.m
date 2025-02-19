function OPD = do_predc(IPD)
%
% Callback function to drive the routine "predc.m"
% Perform Predictive Deconvolution 
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%

OPD = discardprocdata; 
% Trap common errors
if isempty(IPD.d), 
   erh = errordlg('No data to process!', 'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end
% Proceed ...
OPD = IPD; 
[OPD.d, nf, lp] = predc(IPD.d, IPD.dt); 
if isempty(OPD.d),
    disp('PREDICTIVE DECONV > Operation aborted - No data returned!');
    return
end
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab); 
% Update history
iss = size(OPD.history,1); 
text = ['Predictive Deconvolution: Operator length = ' num2str(nf*IPD.dt,4)...
     'ns; Prediction distance = ' num2str(lp*IPD.dt) 'ns.'];
OPD.history(iss+1,1) = cellstr(text); 
return

