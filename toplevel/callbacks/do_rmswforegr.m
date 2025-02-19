function OPD = do_rmswforegr(IPD)
%
% Callback function to drive the routine "rmswforegr.m"
% Remove a sliding window background trace to dipping features
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
[OPD.d, ww] = rmswforegr(IPD.d);   
if isempty(OPD.d),
    disp('RMSWFOREGR > Operation aborted - No data returned!');
    return
end
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab); 
% Update history
iss = size(OPD.history,1); 
fsize = ['[ ' num2str(ww) ' ]'];
text = ['Removed a width ' fsize ' Sliding Window Foreground trace'];
title(text);
OPD.history(iss+1,1) = cellstr(text); 
return
