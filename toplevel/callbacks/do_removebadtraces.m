function OPD = do_removebadtraces(IPD)
%
% Callback function to drive the routine "removebadtraces.m".
% Pick bad traces and replace them with the interpolants computed from
% their near neighbours
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%

OPD = discardprocdata;                           % discard the current OPD
% Trap a quite common mistake 
if isempty(IPD.d),                               
   erh = errordlg('No data to process!', 'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end

% Now Proceed ...
OPD = IPD;                                       % Copy IPD to OPD

% Run removebadtraces.m
[OPD.d, bad] = removebadtraces(IPD.x,IPD.tt2w,IPD.d);

% Check if Operation has been canceled or aborted and issue a message
if isempty(OPD.d),
    disp('REMOVE BAD TRACES > Operation aborted - No O/P data returned!');
    return
end

% Display data
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab);
title([num2str(size(bad,2)) ' bad traces replaced!'])

% Update processing history
iss = size(OPD.history,1); 
text = ['Replaced ' num2str(size(bad,2)) ' bad traces'];
OPD.history(iss+1,1) = cellstr(text); 

return
