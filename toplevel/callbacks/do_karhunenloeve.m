function OPD = do_karhunenloeve(IPD)
%
% Callback function to drive the routine "karhunenloeve.m"
% Smooth data with a low-dimensional approximation using the 
% Karhunen - Loeve transformation 
%
% Copyright (C) 2005, 2008 Andreas Tzanis. All rights reserved.
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
[OPD.d, P] = KarhunenLoeve(IPD.d, IPD.x, IPD.tt2w); 
if isempty(OPD.d),
    disp('KARHUNEN-LOEVE > Operation aborted - No data returned!');
    return
end
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab); 
% Update history
iss = size(OPD.history,1); 
text = ['Applied Karhunen-Loeve Trans. (Used ' num2str(P) ...
    ' largest eigenvectors).'];
OPD.history(iss+1,1) = cellstr(text); 
return
