function OPD = do_fk_filter(IPD)
%
% Callback function to drive the routine "fk_filter.m"
% Perform interactive F-K filtering
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%

OPD = discardprocdata; 
% Trap common errors
if isempty(IPD.d), 
   erh = errordlg('No data to process!', 'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end
if isempty(IPD.dx),     
    erh = errordlg('F-K filtering possible only for equally-spaced data!', ...
        'MATGPR : ERROR'); 
     uiwait(erh);  
     return
end
% Proceed ...
OPD = IPD;
[OPD.d, fktype, vrange] = fk_filter(IPD.d, IPD.dt, IPD.dx);
if isempty(OPD.d),
    disp('F-K FILTER > Operation aborted - No data returned!');
    return
end
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab); 
title('F-K Filtered Data')
% Update history
iss = size(OPD.history,1); 
if fktype == 1,
    text = 'Applied Polygonal Zone-PASS F-K Filter';
end
if fktype == 2,
    text = 'Applied Polygonal Zone-STOP F-K Filter';
end
if fktype == 3,
    text = ['Velocity-Fan PASS F-K Filter,[' ...
        num2str(vrange(1)) ',' num2str(vrange(2)) ']&[' ...
        num2str(vrange(3)) ',' num2str(vrange(4)) ']m/ns'];
end
if fktype == 4,
    text = ['Velocity-Fan STOP F-K Filter,[' ...
        num2str(vrange(1)) ',' num2str(vrange(2)) ']&[' ...
        num2str(vrange(3)) ',' num2str(vrange(4)) ']m/ns'];
end
if fktype == 5,
    text = 'Applied UP-DIP F-K Filter';
end
if fktype == 6,
    text = 'Applied DOWN-DIP F-K Filter';
end
OPD.history(iss+1,1) = cellstr(text); 
return
