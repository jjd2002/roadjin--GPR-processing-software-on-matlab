function OPD = do_k_filter(IPD, filtertype)
%
% Callback function to drive FIR wavenumber filtering 
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

OPD = discardprocdata;  
if isempty(IPD.d),
    erh = errordlg('No data to process!','MATGPR: ERROR');
    uiwait(erh);
    return
end
if isempty(IPD.dx),      
    erh = errordlg(['Wavenumber filtering possible only for ' ... 
        'equally-spaced data!'], 'MATGPR : ERROR'); 
    uiwait(erh); 
    return
end;
OPD = IPD; 
[OPD.d, Ws] = k_filter(IPD.tt2w,IPD.d,IPD.dx,filtertype, ENVAR.Design_Filter); 
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
    text = ['Applied Low-Pass K-Filter at ' ...
        num2str(Ws(1)*fn) ' m^(-1)'];
elseif filtertype == 2,
    text = ['Applied High-Pass K-Filter, at ' ...
        num2str(Ws(1)*fn) ' m^(-1)'];
elseif filtertype == 3,
    text = ['Applied Band-Pass K-Filter, ['...
        num2str(Ws(1)*fn,4) '- ' num2str(Ws(2)*fn,4) '] m^(-1)'];
elseif filtertype == 4,
    text = ['Applied Band-Stop K-Filter, ['...
        num2str(Ws(1)*fn,4) '- ' num2str(Ws(2)*fn,4) '] m^(-1)'];
end
OPD.history(iss+1,1) = cellstr(text); 
title(text);
return
