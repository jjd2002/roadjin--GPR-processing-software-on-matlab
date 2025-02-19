function do_exportdata(IPD,DataFormat)
%
% Callback function to export the IPD structure by driving the routines:
% "writesu.m"    if DataFormat = 'su_format'
% "writesegy.m", if DataFormat = 'segy_format'
% "writedzt.m",  if DataFormat = 'dzt_format'
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%

% Trap common errors ...
if isempty(IPD.d),                               
   erh = errordlg('There is NO DATA to export!', 'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end
if isempty(IPD.dx),
    if strcmpi(DataFormat,'su_format'),
        erh = errordlg(['Only equally-spaced data are exportable ' ...
            'to SU format!'], 'MATGPR: ERROR'); 
    end
    if strcmpi(DataFormat,'segy_format'),
        erh = errordlg(['Only equally-spaced data are exportable ' ...
            'to SEG-Y format!'], 'MATGPR: ERROR'); 
    end
    uiwait(erh); 
    return;
end;
% Proceed ... 
if strcmpi(DataFormat,'su_format'),
    writesu(IPD);                                 % export SU format
end
if strcmpi(DataFormat,'segy_format'),
    writesegy(IPD);                               % export SEG-Y format
end
if strcmpi(DataFormat,'dzt_format'),
    writedzt(IPD);                               % export DZT format
end
return
