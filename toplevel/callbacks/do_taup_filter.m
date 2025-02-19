function OPD = do_taup_filter(IPD)
%
% Callback function to drive the routine "taup_filter.m"
% Perform interactive data modeling and denoising in the TAUP-P domain
%
% Copyright (C) 2008, Andreas Tzanis. All rights reserved.
%

global rdtype

OPD = discardprocdata; 
% Trap common errors
if isempty(IPD.d), 
   erh = errordlg('No data to process!', 'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end
if isempty(IPD.dx),     
    erh = errordlg('TAU-P filtering possible only for equally-spaced data!', ...
        'MATGPR : ERROR'); 
     uiwait(erh);  
     return
end
% Proceed ...
OPD = IPD;
OPD.d = taup_filter(IPD.d, IPD.tt2w, IPD.dt, IPD.x, IPD.dx);
if isempty(OPD.d),
    disp('TAU-P FILTER > Operation aborted - No data returned!');
    return
end
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab); 
title('TAU-P Filtered Data')
% Update history
iss = size(OPD.history,1); 
if rdtype == 11,
    text = 'Applied Zone-PASS Tau-P Filter and kept the MODEL';
end
if rdtype == 12,
    text = 'Applied Zone-PASS Tau-P Filter and kept the RESIDUALS';
end
if rdtype == 21,
    text = 'Applied Zone-STOP Tau-P Filter and kept the MODEL';
end
if rdtype == 22,
    text = 'Applied Zone-STOP Tau-P Filter and kept the RESIDUAS';
end
OPD.history(iss+1,1) = cellstr(text); 
return
