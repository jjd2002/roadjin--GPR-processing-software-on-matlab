function OPD = do_removedc(IPD)
%
% Callback function to drive the routine "removedc.m"
% Remove DC from each trace individually
%
% Copyright (C) 2006, Andreas Tzanis. All rights reserved.
%

OPD = discardprocdata;                           % discard the current OPD
% Trap common errors ...
if isempty(IPD.d),                               
   erh = errordlg('请先导入数据!', 'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end
% Now Proceed ...
OPD = IPD;                                       % Copy IPD to OPD
% Run removedc.m
OPD.d = removedc(IPD.d);

% Check if Operation has been canceled or aborted and issue a message
if isempty(OPD.d),
    disp('REMOVEDC > Operation aborted - No O/P data returned!');
    return
end

% Display data
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab);
title('Removed DC component from individual traces')

% Update processing history
iss = size(OPD.history,1); 
text = 'Removed DC component';
OPD.history(iss+1,1) = cellstr(text); 

return
