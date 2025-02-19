function OPD = do_gainagc(IPD)
%
% Callback function to drive the routine "gainagc.m"
% Standard Automatic Gain Control
%
% Copyright (C) 2005, 2010 Andreas Tzanis. All rights reserved.
%

OPD = discardprocdata; 
% Trap common errors
if isempty(IPD.d), 
   erh = errordlg('请先导入数据!', 'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end
% Proceed ...
OPD = IPD;
% Get AGC window length in ns
ask    = inputdlg('给定时窗 (ns)','GAINAGC : REQUEST',1);
if isempty(ask),
    OPD = discardprocdata;
    disp('GAIN AGC > Operation aborted by user!');
    return
end
wagc   = str2num(ask{1});                      % agc window in ns
iwagc = floor(wagc/IPD.dt);                    % agc window in samples
if iwagc < 1, 
    erh = errordlg([ 'AGC 时窗要大于1'],'GAINAGC : ERROR');
    uiwait(erh)
    OPD = discardprocdata;
    return
end
if iwagc > IPD.ns/2, 
    erh = errordlg([ 'AGC 时窗 = ' num2str(wagc) '太长了! '],...
        'GAINAGC : ERROR');
    uiwait(erh)
    OPD = discardprocdata;
    return
end
% Apply AGC
OPD.d = gainagc(IPD.d,IPD.dt, wagc); 
if isempty(OPD.d),
    disp('GAIN AGC > Operation aborted - No data returned!');
    return
end
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab); 
% Update history
iss = size(OPD.history,1); 
text = 'Standard AGC';
title(text)
OPD.history(iss+1,1) = cellstr(['Applied ' text]); 
return
