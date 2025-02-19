function OPD = do_editscanaxis(IPD);
%
% Callback function to drive the routine "editscanaxis.m"
% Edit the scan axis (handle groups of traces)
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
[OPD.d,OPD.x,OPD.ntr,OPD.markertr,OPD.xyz] =  ...
         editscanaxis(IPD.d,IPD.x,IPD.dx,IPD.markertr,IPD.xyz); 
if isempty(OPD.d),
    disp('EDIT SCAN AXIS > Operation aborted - No data returned!');
    return
end
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab); 
% Update history
iss = size(OPD.history,1); 
text = ['Trimmed Scan Axis; cleared ' num2str(IPD.ntr-OPD.ntr) ...
          ' traces.'];
OPD.history(iss+1,1) = cellstr(text); 
return
