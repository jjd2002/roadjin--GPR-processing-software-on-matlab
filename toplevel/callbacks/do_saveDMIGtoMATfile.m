function do_saveDMIGtoMATfile(OPD)
%
% Callback function to save in a MAT file, a depth migrated section (that
% started life as an OPD structure). 
%
% Copyright (C) 2005,2008 Andreas Tzanis. All rights reserved.
%

global ENVAR
% Trap a common error
if isempty(OPD.d), 
    erh = errordlg('No data to save in the OPD structure!', ...
        'MATGPR: ERROR'); 
    uiwait(erh);  
    return; 
end; 
if ~isempty(findstr(OPD.zlab,'Depth')) && ~isempty(OPD.dz),
    getts  = findstr(OPD.fname,'.');
    DMIGfname = [OPD.fname(1:getts-1) '_dmig.' ENVAR.binary_file_type];
    if strcmp(ENVAR.binary_file_type,'mgp'),
        Ofiles = {'*.mgp;', 'MGP-files (*.mgp)';
            '*.mat','MAT-files (*.mat)'; ...
            '*.mgp', 'MGP-files (*.mgp)'};
    elseif strcmp(ENVAR.binary_file_type,'mat'),
        Ofiles = {'*.mat;', 'MAT-files (*.mat)';
            '*.mgp', 'MGP-files (*.mgp)'; ...
            '*.mat','MAT-files (*.mat)'};
    end
    [DMIGfname, OPD.pname]= uiputfile( Ofiles, ...
        'Save depth migrated data as ...',[ENVAR.currentworkdir DMIGfname]);
    if DMIGfname == 0;     % cancel was pressed
        return;
    end;
    OPD.fname = DMIGfname;
    IPD            = OPD;
    IPD.dt         = [];
    IPD.tt2w       = [];
    IPD.TimesSaved = 0;
    if strcmp(ENVAR.binary_file_type, 'mat'),
        save([IPD.pname DMIGfname],'IPD',ENVAR.matfileformat);
    elseif strcmp(ENVAR.binary_file_type, 'mgp'),
        save2MGP(IPD,[IPD.pname DMIGfname]);
    end
else
    erh = errordlg('I do not detect any depth migrated data!', ...
        'MATGPR: ERROR');
    uiwait(erh);
end
return
