function IPD = do_importdata(DataFormat)
%
% Callback function to import GPR data by driving the routines:
% "readdzt.m"            if DataFormat = 'dzt_format'
% "readrd3.m"            if DataFormat = 'rd3_format'
% "readdt1.m"            if DataFormat = 'dt1_format'
% "readZONDsegydata.m"   if DataFORMAT = 'zond_segy'
% "readsu.m"             if DataFormat = 'su_format'
% "readsegy.m",          if DataFormat = 'segy_format'
% "readhcd.m",          if DataFormat = 'hcd_format'
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%

global ENVAR
% Discard the current OPD
OPD = discardprocdata; 

% Proceed ... 
if strcmpi(DataFormat,'dzt_format'),
    IPD = readdzt;                                % get GSSI DZT data
end
if strcmpi(DataFormat,'rd3_format'),
    IPD = readrd3;                                % get Mala RD3 data
end
if strcmpi(DataFormat,'dt1_format'), 
    IPD = readdt1;                                % get PULSE EKKO data
end
if strcmpi(DataFormat,'zond_segy'),
    IPD = readZONDsegydata;                       % get ZOND data
end
if strcmpi(DataFormat,'su_format'),
    IPD = readsu;                                 % get Seismic Unix format
end
if strcmpi(DataFormat,'segy_format'),
    IPD = readsegy;                               % get SEGY format
end
if strcmpi(DataFormat,'hcd_format'),
    IPD = readhcd;                               % get HCD format
end
% Check if operation canceled or aborted
if isempty(IPD.d),
   erh = warndlg(['请选择数据!'], 'MATGPR: WARNING'); 
   uiwait(erh);  
   return
end

% When importing raw data, always save to the preferred binary format
if strcmp(ENVAR.binary_file_type, 'mat'),
    IPD.fname = [IPD.fname '.mat'];
    save([IPD.pname IPD.fname],'IPD',ENVAR.matfileformat);
elseif strcmp(ENVAR.binary_file_type, 'mgp'),
    IPD.fname = [IPD.fname '.mgp'];
    save2MGP(IPD,[IPD.pname IPD.fname]);
end

% Erase the UNDO/RESTORE buffer of the previous session, if any
set(findobj('tag','fi0'),'userdata',[]);
set(findobj('tag','undoxx'),'Enable','off');
set(findobj('tag','undoundo'),'Enable','off');

% Display IPD information on the MATGPR window
showinfo(IPD);

% Display the data
viewdata(IPD.x,IPD.tt2w,IPD.d,'indata',IPD.xlab,IPD.zlab); 

% Reset the current work directory
ENVAR.currentworkdir = IPD.pname;

return
