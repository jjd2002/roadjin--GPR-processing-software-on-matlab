function OPD = do_equalize(IPD)
%
% Callback function to drive the routine "equalize.m"
% Equalize traces to have a sum of absolute values equal to some base value
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%

OPD = discardprocdata;                           % discard the current OPD
% Trap common errors ...
if isempty(IPD.d),                               
   erh = errordlg('No data to process!', 'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end

% Now Proceed ...
OPD = IPD;                                       % Copy IPD to OPD

% Run SOMETHING.m
OPD.d = equalize( IPD.d );

% Check if Operation has been canceled or aborted and issue a message
if isempty(OPD.d),
    disp('EQUALIZE > Operation aborted - No O/P data returned!');
    return
end

% Display data
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab);
title('GPR section with EQUALIZED traces')

% Update processing history
iss = size(OPD.history,1); 
text = 'Equalized all traces';
OPD.history(iss+1,1) = cellstr(text); 

return
