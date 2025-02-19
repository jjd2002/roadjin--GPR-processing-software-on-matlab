function IPD = concatenate(action)
%
%CONCATENATE : Join or concatenate two or more GPR radar data sets. 
%==>   Unequally spaced data (no survey wheels and prior to marker
%      interpolation) - IPD.dx is empty: This type of data can be handled
%      but all necessary trace editing and trimming operations must be done
%      before calling CONCATENATE. Marker trace information is preserved
%      and the concatenated section can be subsequently interpolated to
%      equal spacing.  
%==>   Equally spaced data (with survey wheels or after marker
%      interpolation (IPD.dx assigned): All data sets must be in MATGPR
%      (.mat) format and have the same number of samples per trace
%      (IPD.ns), the same sampling rate (IPD.dt) and the same trace spacing
%      (IPD.dx). Marker trace information is preserved in the resulting
%      section, but it is adviseable that trace coordinates (IPD.xyz) be
%      assigned before CONCATENATE is called, otherwise IPD.xyz will be
%      empty.  
%==>   For GSSI (.DZT)data, header gain manipulations must have been done
%      before CONCATENATE is called, because IPD.DZThdgain is not
%      preserved.  
% 
% Usage : concatenate('initialize') to set up the GUI for assembling the
%         data files to concatenate. Afterwards, the routine recurses on
%         itself (IPD = concatenate('process') to perform the
%         concatenation. 
%
%Author : Andreas Tzanis,
%         Department of Geophysics, 
%         University of Athens
%         atzanis@geol.uoa.gr
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%

global files2concat     % assembly of sections (files) to concatenate

% Startup the GUI to assemble profiles for concatenation
if strcmp(action,'initialize'),
    
% Create GUI interface
    figure('Name','Concatenate Sections', 'Tag', 'concat_gui', ...
        'NumberTitle','off', 'Position',[200 100 550 500 ], ...
        'MenuBar','none', 'Color',[0 0.5 0.5] );
    axis off
% Get filenames to concatenate
    text('Color',[1 1 1 ],  'Units','normalized', ...
        'Position', [0.0 0.93 0 ], 'FontSize',10 , ...
        'FontWeight','bold', 'String','Path / File:');
    uicontrol('Style','edit', 'Tag', 'filebox', ... 
        'Units','normalized', 'Position',[.05 .80 .75 .045], ...
        'BackGroundColor','w', 'String',' ', ... 
        'Callback', @filebx);
    uicontrol('Style','Pushbutton','tag','browsebutn', ...
        'Units','normalized', 'Position',[.82 .80 .15 .045 ], ...
        'BackGroundColor','w','String','Browse', 'Callback', @brse);
% Store paths and file names in a list box
    uicontrol('Style', 'listbox','Tag','lsbox', ...
        'Units','normalized', 'Position',[0.05 0.2 0.75 0.55], ...
        'string', cell(0,0));
% Remove files from listbox
    uicontrol('Style','Pushbutton','tag','rmlsbutn', ...
        'Units','normalized', 'Position',[.82 .705 .15 .045 ], ...
        'BackGroundColor','w','String','Remove from List', ...
        'Callback', @rmfilefromlist);
% Proceed ...
    uicontrol('Style','Pushbutton','tag','gobutn', ...
        'Units','normalized', 'Position',[.82 .475 .15 .045 ], ...
        'BackGroundColor','w','String','GO', ...
        'Callback', 'IPD=concatenate(''process''); ');
% Cancel ...
    uicontrol('Style','Pushbutton','Units','normalized', ... 
        'Position',[.82 .2 .15 .045 ], 'BackGroundColor','w', ...
        'String','Cancel', 'Callback', ...
        'delete(findobj(''tag'',''concat_gui'')); return; ');
    return
end

% File assembly collected - now concatenate
if strcmp(action,'process'),
    
    files2concat = get(findobj('tag','lsbox'),'string');
    lf2c = length(files2concat);
    sps = zeros(lf2c,1);
    dtt = zeros(lf2c,1);
    spm = zeros(lf2c,1);
    tro = zeros(lf2c,1);
    xyz = zeros(lf2c,1);
% First pass - check consistency of parameter sets 
    for i=1:lf2c,
        if ~isempty(findstr(files2concat{i},'.mat')),
            IPD = importdata(files2concat{i});
        end
        if ~isempty(findstr(files2concat{i},'.mgp')),
            IPD = readMGP(files2concat{i});
        end

        sps(i) = IPD.ns;
        dtt(i) = IPD.dt;
        if isempty(IPD.dx),
            spm(i) = 0;
        else
            spm(i) = IPD.dx;
        end
        tro(i) = IPD.TxRx;
        if ~isempty(IPD.xyz.Tx),
            xyz(i) = 1;
        end
        if i==1,
            point1 = IPD.x(1);
        end
    end
% Error trapping
    if ~isempty(find(diff(dtt)~=0)),
        erh = errordlg(['At least 2 sections have different sampling '...
            'rates! Aborting!'],'CONCATPROFILES: ERROR');
        uiwait(erh);
        IPD = initdatastr;
        delete(findobj('tag','concat_gui')); 
        return
    end
    if ~isempty(find(diff(sps)~=0)),
        erh = errordlg(['At least 2 sections have different '...
            'traveltime vectors! Aborting!'],'CONCATPROFILES: ERROR');
        uiwait(erh);
        IPD = initdatstr;  
        delete(findobj('tag','concat_gui')); 
        return
    end
    if ~isempty(find(diff(spm)~=0)),
        erh = errordlg(['At least 2 sections have different trace '...
            'spacings! Aborting!'],'CONCATPROFILES: ERROR');
        uiwait(erh);
        IPD = initdatastr;  
        delete(findobj('tag','concat_gui')); 
        return
    end
    if ~isempty(find(diff(tro)~=0)),
        erh = warndlg(['At least 2 sections have different antenna '...
            'separations! Proceed at your own risk!'], ...
            'CONCATPROFILES: WARNING');
        uiwait(erh);
    end
    xyzOK = 1;
    if ~isempty(find(diff(xyz)~=0)),
        erh = warndlg(['Some sections have XYZ data and some have not. '...
            'XYZ data will be excluded from the concatenated section!'], ...
            'CONCATPROFILES: WARNING');
        uiwait(erh);
        xyzOK = 0;
    end
% ALL is OK - concatenate sections
    OPD = initdatastr;
    for i=1:lf2c,
        if ~isempty(findstr(files2concat{i},'.mat')),
            IPD = importdata(files2concat{i});
        end
        if ~isempty(findstr(files2concat{i},'.mgp')),
            IPD = readMGP(files2concat{i});
        end

        OPD.d   = [OPD.d IPD.d];
        ntrtmp = size(OPD.d,2);
        if ~isempty(IPD.markertr),
            IPD.markertr(:,1) = IPD.markertr(:,1) + ntrtmp - IPD.ntr;
            OPD.markertr = [OPD.markertr; IPD.markertr];
        end
        if xyzOK,
            OPD.xyz.Tx = [OPD.xyz.Tx; IPD.xyz.Tx];
            OPD.xyz.Rx = [OPD.xyz.Rx; IPD.xyz.Rx];
        end
    end
    OPD.dt    = IPD.dt;
    OPD.ns    = IPD.ns;
    OPD.tt2w  = IPD.tt2w;
    OPD.ntr   = size(OPD.d,2);
    if isempty(IPD.dx),
        OPD.dx   = [];
        OPD.x    = 1 : 1 : OPD.ntr;
        OPD.xlab = 'Scan Line (# Traces)';
    else
        OPD.dx   = IPD.dx;
        OPD.x    = point1 : OPD.dx : (OPD.ntr-1)*OPD.dx;
        OPD.xlab = 'Scan Line (meters)';
    end
% Confirm TxRx offset
    answer  = inputdlg('Confirm Tx-Rx offset:','CONCATPROFILES: Request',...
        1,{num2str(IPD.TxRx)});
    if isempty(answer),
        OPD.TxRx = 0; 
    else
        OPD.TxRx = str2num(answer{1});       
    end
% Confirm Antenna type
    answer  = inputdlg('Confirm antenna model:','CONCATPROFILES: Request',...
        1,{IPD.Antenna});
    if isempty(answer),
        OPD.Antenna  = 'Unknown';
    else
        OPD.Antenna  = answer{1};       
    end
% Remaining parameters
    OPD.origin    = [IPD.origin  ' (Concatenated Sections - see comments)'];
    OPD.comments  = strvcat(['This data set was created by ' ...
        'concatenating the GPR sections:'], str2mat(files2concat));
    OPD.sigpos    = 0;
    OPD.DZThdgain = [];     % Gain operations must have been done by now 
    OPD.zlab      = 'Traveltime (ns)';
    IPD = OPD;              % assign output data structure
% Wrap up and exit
    delete(findobj('tag','concat_gui')); 
    return
end

function filebx(filebox,eventdata)
global files2concat
xfname = get(findobj('tag','filebox'),'string'); 
if exist(xfname)==2, 
    files2concat = get(findobj('tag','lsbox'),'string'); 
    nfiles = length(files2concat); 
    files2concat(nfiles+1) = cellstr(xfname);  
    set(findobj('tag','lsbox'),'string',files2concat); 
else
    erw = warndlg('No such file!','Warning'); 
    uiwait(erw); 
end
return

function brse(browsebtn,eventdata)
global ENVAR files2concat
[xfname, xpname]= uigetfile( ...
    {'*.mat; *.mgp','MATGPR Files (*.mat, *.mgp)';
     '*.mat','MAT-files (*.mat)'; ...
     '*.mgp', 'MGP-files (*.mgp)'}, ...
     'Select Files',ENVAR.currentworkdir,'multiselect','on'); 
if isempty(xfname),                                   % Operation cancelled
    return
end
files2concat = get(findobj('tag','lsbox'),'string');
nfiles = length(files2concat);
if iscellstr(xfname)                              % Multiple files selected
    for i=1:max(size(xfname))
        files2concat(nfiles+i) = cellstr([xpname xfname{i}]);
    end
    set(findobj('tag','lsbox'),'string',files2concat);
elseif ischar(xfname)                                % Single file selected 
    set(findobj('tag','filebox'),'string',[xpname xfname]); 
    files2concat(nfiles+1) = cellstr([xpname xfname]);  
    set(findobj('tag','lsbox'),'string',files2concat); 
end
return

function rmfilefromlist(rmlsbutn,eventdata)
global files2concat
nf2c = get(findobj('tag','lsbox'),'value');  
files2concat = get(findobj('tag','lsbox'),'string');  
lf2c = length(files2concat);  
if nf2c == 1,  
    files2concat = files2concat(2:lf2c);  
elseif nf2c == lf2c,  
    files2concat = files2concat(1:lf2c-1);  
else
    files2concat = [files2concat(1:nf2c-1); ...
        files2concat(nf2c+1:length(files2concat))];  
end  
set(findobj('tag','lsbox'),'value', 1,'string',files2concat);  
return

