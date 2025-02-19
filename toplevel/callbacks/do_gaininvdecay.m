function OPD = do_gaininvdecay(IPD)
%
% Callback function to drive the Inverse Decay gain function 
%
% Copyright (C) 2009, Andreas Tzanis. All rights reserved.
%

OPD = discardprocdata; 
if isempty(IPD.d), 
   erh = errordlg('No data to process!', 'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end
OPD = IPD;
OPD.d = gaininvdecay(IPD.d, IPD.tt2w); 
if isempty(OPD.d), 
    disp('GAININVDECAY > Operation aborted - No O/P data returned!');
    return
end
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab); 
iss = size(OPD.history,1); 
text = 'Applied Gain by Inverse Amplitude Decay';
OPD.history(iss+1,1) = cellstr(text); 
return
