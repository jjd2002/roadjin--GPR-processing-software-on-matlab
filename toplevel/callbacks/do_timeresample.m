function OPD = do_timeresample(IPD)
%
% Callback function to drive routine "resample1.m"
% Resample the time axis at a higher or lower rate
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%

OPD = discardprocdata; 
% Trap a most common mistake 
if isempty(IPD.d),                               
   erh = errordlg('No data to process!', 'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end
% Proceed ..
OPD = IPD;
[OPD.d, OPD.dt, OPD.ns] = timeresample(IPD.d, IPD.dt);  

% Check if operation has been canceled or aborted and issue a message
if isempty(OPD.d),
    disp('TIME RESAMPLE > Operation aborted - No O/P data returned!');
    return
end

% Update O/P data traveltime vector
OPD.tt2w = 0 : OPD.dt : (OPD.ns -1)*OPD.dt; 
% If DZT data, update header gain information
if ~isempty(IPD.DZThdgain), 
    OPD.DZThdgain = interp1(IPD.tt2w,IPD.DZThdgain,OPD.tt2w);
end

% Display processed data
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab); 

% Display title and update processing history
text = ['Resampled Time Axis, from ' num2str(IPD.ns) ...
       ' to ' num2str(OPD.ns) ' samples.'];
title(text);
iss = size(OPD.history,1); 
OPD.history(iss+1,1) = cellstr(text); 
return
