function OPD = do_statics(IPD)
%
% Callback function to drive the routine "staticsgui.m"
% Perform topographic corrections to a GPR section                  
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%

OPD = discardprocdata;                           % discard the current OPD
% Trap common errors ...
% A most common mistake 
if isempty(IPD.d),                               
   erh = errordlg('No data to process!', 'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end
% Another common mistake - delete if not applicable
if isempty(IPD.xyz.Tx),                              
    erh = errordlg(['Topographic information is MISSING. Please check', ...
        ' and try again!'], 'MATGPR : ERROR');
    uiwait(erh); 
    return;
end;
% Now Proceed ...
OPD = IPD;                                       % Copy IPD to OPD
% Run staticsgui.m
OPD.d = staticsgui(IPD.d, IPD.dt, IPD.xyz.Tx, IPD.xyz.Rx); 
% Check if Operation has been canceled or aborted and issue a message
if isempty(OPD.d),
    disp('DO_STATICS > Operation aborted - No O/P data returned!');
    return
end
% Display data
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab);
title('Applied Topographic Correction')
% Update processing history
iss = size(OPD.history,1); 
text = 'Applied Topographic Corrections';
OPD.history(iss+1,1) = cellstr(text); 

return


