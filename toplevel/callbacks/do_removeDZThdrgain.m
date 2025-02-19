function OPD = do_removeDZThdrgain(IPD)
%
% Callback function to drive DZT header range-gain removal 
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%

OPD = discardprocdata; 
if ~isempty(IPD.DZThdgain),     
    OPD = IPD; 
    OPD.d = gainrmdzthdr(IPD.d, IPD.DZThdgain); 
    if ~isempty(OPD.d), 
        viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab); 
        iss = size(OPD.history,1); 
        OPD.history(iss+1,1) = cellstr('Removed DZT header gain');
    end; 
else 
    erh = errordlg('Are you sure this is GSSI (DZT) data ?', ...
            'MATGPR: ERROR'); 
    uiwait(erh);   
    return; 
end;   
return
