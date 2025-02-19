function OPD = do_dewow(IPD)
%
% Callback function to drive the routine "dewow.m"
% Apply a dewow filter
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
OPD.d = dewow(IPD.d);   
if isempty(OPD.d),
    disp('DEWOW > Operation aborted - No data returned!');
    return
end
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab); 
% Update history
iss = size(OPD.history,1); 
text = 'Applied Dewow Filter';
OPD.history(iss+1,1) = cellstr(text); 
return
