function IPD = do_readdzt()
%
% Callback function to drive the routine "readdzt.m"
% Import raw GSSI (DZT) data
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%

OPD = discardprocdata;                         % discard the current OPD
IPD = readdzt;                                 % get DZT data
if ~isempty(IPD.d),                            % test if operation canceled
    IPD.fname = [IPD.fname '.mat'];
    save([IPD.pname IPD.fname],'IPD');         % save to MAT file
    showinfo(IPD);                             % display information
    viewdata(IPD.x,IPD.tt2w,IPD.d,'indata',IPD.xlab,IPD.zlab); 
end
return
