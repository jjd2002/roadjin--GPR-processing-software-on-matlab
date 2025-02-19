function OPD = do_gaintpow(IPD)
%
% Callback function to drive the scale*t^pow gain function 
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%

OPD = discardprocdata; 
if isempty(IPD.d), 
   erh = errordlg('请先导入数据!', 'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end
OPD = IPD;
[OPD.d, scale, pow] = gaintpow(IPD.d, IPD.tt2w); 
if isempty(OPD.d), 
    disp('GAINTPOW > Operation aborted - No O/P data returned!');
    return
end
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab); 
iss = size(OPD.history,1); 
text = [ 'Applied gain g(t)= ' num2str(scale,4) '*t^' num2str(pow,3)];
OPD.history(iss+1,1) = cellstr(text); 
return
