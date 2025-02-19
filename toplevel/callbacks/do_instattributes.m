function OPD = do_instattributes(IPD,attribute);
%
% Callback function to compute and view instantaneous attributes of the
% input data.
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%

OPD = discardprocdata; 
OPD = IPD; 
OPD.d = instattributes(IPD.d, IPD.dt, attribute);
if ~isempty(OPD.d); 
    viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab);
    iss = size(OPD.history,1); 
    if strcmp(lower(attribute),'amplitude'),
        text = [' Instantaneous ' attribute ];
    elseif strcmp(lower(attribute),'atan'),
        text = [' Instantaneous  Phase [-90, 90]' ];
    elseif strcmp(lower(attribute),'atan2'),
        text = [' Instantaneous  Phase [-180, 180]' ];
    elseif strcmp(lower(attribute),'unwraped'),
        text = [' Instantaneous  Phase, Continuous (unwraped)' ];
    elseif strcmp(lower(attribute),'ifreq'),
        text = [' Instantaneous  Frequency' ];
    end
    OPD.history(iss+1,1) = cellstr(text);
% Give a title to the "outdata" figure
    title(text,'fontsize',12);
end 
return
