function OPD = do_f_filter(IPD, filtertype)
%
% Callback function to drive FIR frequency filtering 
%
% filtertype : flag determining the type of the filter
%               =1  Low-pass filter
%               =2  High-pass filter
%               =3  Band-pass filter
%               =4  Band-stop filter
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%

global ENVAR
OPD = discardprocdata;                           % discard the current OPD
% Trap common errors ...
if isempty(IPD.d),                               
   erh = errordlg('No data to process!', 'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end
% Now Proceed ...
OPD = IPD;                                       % Copy IPD to OPD
% Run f_filter
[OPD.d, Ws] = f_filter(IPD.x,IPD.d,IPD.dt,filtertype, ...
                      ENVAR.Design_Filter); 
% No result for some reason ...
if isempty(OPD.d),
    disp('F_FILTER > Operation aborted - No O/P data returned!');
    return
end
% Display results
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab);
% Update processing history
iss = size(OPD.history,1); 
fn = 1/(2*IPD.dt);
if filtertype == 1,
 text = ['Applied Low-Pass F-Filter at ' ...
     num2str(Ws(1)*fn) ' GHz'];
elseif filtertype == 2,
 text = ['Applied High-Pass F-Filter, at ' ...
     num2str(Ws(1)*fn) ' GHz'];
elseif filtertype == 3,
 text = ['Applied Band-Pass F-Filter, ['...
     num2str(Ws(1)*fn,4) '- ' num2str(Ws(2)*fn,4) '] GHz'];
elseif filtertype == 4,
 text = ['Applied Band-Stop F-Filter, ['...
     num2str(Ws(1)*fn,4) '- ' num2str(Ws(2)*fn,4) '] GHz'];
end
OPD.history(iss+1,1) = cellstr(text); 
title(text);
return
