function OPD = do_gaingagc(IPD);
%
% Callback function to drive the Gaussian tapered AGC program 
%
% Copyright (C) 2005, 2010, Andreas Tzanis. All rights reserved.
%

OPD = discardprocdata; 

if isempty(IPD.d), 
   erh = errordlg('请先导入数据!', 'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end
% Proceed ...
OPD = IPD;
% Get AGC parameters
ask       = inputdlg({'定义噪音水平（例如，5.0e-7)' ...
    '给定时窗(ns)'},...
    'Gaussian-tapered AGC ',1);
if isempty(ask), 
    OPD = discardprocdata;
    disp('GAUSSIAN TAPERED AGC > Operation aborted by user!');
    return; 
end;
% Noise Level (Gaussian Window Width)
NoiseLevel = str2num(ask{1});
if ~NoiseLevel,
    NoiseLevel = 5.0e-7;
end
% AGC window size
wagc      = str2num(ask{2});
if ~wagc,
    erh = errordlg('Wrong data specification','GT-AGC : ERROR');
    uiwait(erh)
    OPD = discardprocdata;
    return
end
iwagc  = floor(wagc/IPD.dt);                       % AGC window in samples
if iwagc < 1, 
    erh = errordlg([ 'AGC 时窗要大于1'],'GAINGAGC : ERROR');
    uiwait(erh)
    OPD = discardprocdata;
    return
end
if iwagc > IPD.ns/2, 
    erh = errordlg([ 'AGC window = ' num2str(iwagc) '太长了! '],...
        'GAINGAGC : ERROR');
    uiwait(erh)
    OPD = discardprocdata;
    return
end
% Apply the GT-AGC
OPD.d = gaingagc(IPD.d,IPD.dt, NoiseLevel, wagc); 
if isempty(OPD.d),
    disp('GAUSSIAN-TAPERED AGC > Operation aborted - No data returned!');
    return
end
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab);
iss = size(OPD.history,1); 
text = 'Gaussian-tapered AGC';
title(text);
OPD.history(iss+1,1) = cellstr(['Applied ' text]); 
return
