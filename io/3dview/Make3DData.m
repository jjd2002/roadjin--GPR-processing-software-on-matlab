function D3D = Make3DData(action)
%
% MAKE3DDATA: Create a 3-D volume of GPR data sets for% three-dimensional
%             presentation.    
%
% 1. Presupposes the following Co-ordinate system:
%          ^
%          |
%          |...........................................>
%          |    Direction of survey (scan) lines 
%   Y-axis |...........................................>
%          |
%          |...........................................>
%          |
%          |____________________________________________>
%                  Scan Axis (X-axis)
%
%    ==> The survey is conducted along parallel lines (profiles)
%    ==> The X-axis is the longitudinal direction of the profiles 
%    ==> The Y-axis os the meridional direction, i.e. perpendicular to the
%        profiles.
%
% 2. Other pre-requisites:
%    ==> The only acceptable data format are the MATGPR-specific formats
%        (MGP-files or MAT-files)
%    ==> It is imperative that the X- and Y- coordinates of all traces in
%        all profiles have been assigned a priori - See "Make XYZ" item of
%        the "Basic Handling" menu in the MATGPR documentation.
% 
% 3. Other comments
%    ==> The vertical axis can be either time or depth.
%    ==> The profiles need not have the same lengths, or trace spacings, or
%        traveltime vectors, or sampling rates, or depth vectors, or depth
%        spacings.
%    ==> The profiles need not be equaly spaced along the Y-axis
%
% Usage : Make3DData('initialize') to set up the GUI to assemble the
%         data files from which to create the 3-D volume. Afterwards, the
%         routine recurses on itself ([D3D = Make3DData('process')
%         to create the 3-D volume. 
%
% Output: D3D - Data structure with fields
%         D3D.x, D3D.y, D3D.z : the coordinates of the 3-D data voxels.
%         D3D.d               : the 3-D data volume.
%
%Author : Andreas Tzanis,
%         Department of Geophysics, 
%         University of Athens
%         atzanis@geol.uoa.gr
%
% Copyright (C) 2008, Andreas Tzanis. All rights reserved.
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
    figure('Name','Create 3-D Data', 'Tag', 'make3d_gui', ...
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
        'Callback', 'D3D = Make3DData(''process''); ');
% Cancel ...
    uicontrol('Style','Pushbutton','Units','normalized', ... 
        'Position',[.82 .2 .15 .045 ], 'BackGroundColor','w', ...
        'String','Cancel', 'Callback', ...
        'delete(findobj(''tag'',''make3d_gui'')); return; ');
    return
end

% File assembly collected - now concatenate
if strcmp(action,'process'),
    
    files2concat = get(findobj('tag','lsbox'),'string');
    lf2c   = length(files2concat);
    dat    = cell(1,lf2c);
    ns     = zeros(1,lf2c);
    dtt    = zeros(1,lf2c);
    twt    = cell(1,lf2c);
    dzz    = zeros(1,lf2c);
    zz     = cell(1,lf2c);
    scanax = cell(1,lf2c);
    dx     = zeros(1,lf2c);
    ntrs   = zeros(1,lf2c);
    tro    = zeros(1,lf2c);
    xdata  = cell(1,lf2c);
    ydata  = cell(1,lf2c);
    xyz    = zeros(1,lf2c);
% First pass - make sure that only data with assigned X- and Y- trace
% coordinates may be used!
    for i = 1:lf2c,
        if ~isempty(findstr(files2concat{i},'.mat')),
            IPD = importdata(files2concat{i});
        end
        if ~isempty(findstr(files2concat{i},'.mgp')),
            IPD = readMGP(files2concat{i});
        end
        if ~isempty(IPD.xyz.Tx),
            xyz(i) = 1;
        end
    end

% Second pass - Load data and check consistency of parameter sets 
    i = 0;
    for ii = 1:lf2c,
        if xyz(ii)
            if ~isempty(findstr(files2concat{ii},'.mat')),
                IPD = importdata(files2concat{ii});
            end
            if ~isempty(findstr(files2concat{ii},'.mgp')),
                IPD = readMGP(files2concat{ii});
            end
            i = i+1;
            dat{1,i}    = IPD.d;
            ns(i)      = IPD.ns;
            if ~isempty(IPD.dt)
                twt{1,i}    = IPD.tt2w;
                dtt(i)      = IPD.dt;
            end
            if ~isempty(IPD.dz)
                zz{1,i}    = IPD.z;
                dzz(i)     = IPD.dz;
            end
            if isempty(IPD.dx),
                dx(i) = 0;
            else
                dx(i) = IPD.dx;
            end
            scanax{1,i} = IPD.x;
            ntrs(i)    = IPD.ntr;
            tro(i)     = IPD.TxRx;
            if ~isempty(IPD.xyz.Tx),
                xyz(i) = 1;
                xdata{1,i} = IPD.xyz.Tx(:,1);
                ydata{1,i} = IPD.xyz.Tx(:,2);
            end
        end
    end
    lf2c = i;
 % Check consistency of sampling rates and resample if necessary
 if ~isempty(find(dtt)), %#ok<EFIND>
     if ~isempty(find(diff(dtt)~=0)), %#ok<EFIND>
         disp('At least 2 sections have different sampling rates! ');
         disp('Resampling - new rate will be the median of sampling rates')
         newdt = median(dtt);
         for i=1:lf2c
             oldtwt = twt{1,i};
             mint = oldtwt(1,1);
             maxt = oldtwt(1,ns(i));
             newtwt = (mint:newdt:maxt)';
             newns  = length(newtwt);
             cdat   = dat{1,i};
             [m,n]  = size(cdat);
             cdati  = zeros(newns,n);
             disp(['Resampling section ' num2str(i) '  of  ' num2str(lf2c)]);
             for jj = 1:n
                 cdati(:,jj) = interp1(oldtwt,cdat(:,jj),newtwt);
             end
             dat{1,i} = cdati;
             twt{1,i} = newtwt;
             ns(i) = newns;
             dtt(i) = newdt;
         end
     end
     % Check consistency of traveltime vector lengths and trim or pad with
     % zeros 
     if ~isempty(find(diff(ns)~=0)), %#ok<EFIND>
         oops = menu({'At least two sections have different', 'traveltime vectors!'}, ...
             'Zero-pad all to longest traveltime','Trim all to shortest traveltime',...
             'Abort');
         if oops >= 3 || oops <= 0,
             disp('Operation Aborted!')
             return
         elseif oops == 1,
             newns = max(ns);
             disp('Padding with zeros')
         elseif oops == 2
             newns = min(ns);
         end
         for i=1:lf2c
             cdat    = dat{1,i};
             [m, n]  = size(cdat);
             cdati   = zeros(newns, n);
             if oops == 1,
                 cdati(1:m, 1:n) = cdat;
             elseif oops == 2
                 cdati(1:newns, 1:n) = cdat(1:newns, 1:n);
             end
             dat{1,i}    = cdati;
             ns(i) = newns;
             tt2w = twt{1,i};
             tt2wi = zeros(1,newns);
             if oops == 1,
                 tt2wi(1:m) = tt2w;
             elseif oops == 2,
                 tt2wi(1:newns) = tt2w(1:newns);
             end
             twt{1,i} = tt2wi;
             disp(['Padding section ' num2str(i) '  of  ' num2str(lf2c)]);
         end
     end
 end
 % Check consistency of depth spacings and resample if necessary
 if ~isempty(find(dzz)), %#ok<EFIND>
     if ~isempty(find(diff(dzz)~=0)), %#ok<EFIND>
         disp('At least 2 sections have different depth spacings! ');
         disp('Resampling - new spacing will be the median of depth spacings')
         newdz = median(dzz);
         for i=1:lf2c
             oldzz = zz{1,i};
             minz = oldzz(1,1);
             maxz = oldzz(1,ns(i));
             newzz = (minz:newdz:maxz)';
             newns  = length(newzz);
             cdat   = dat{1,i};
             [m,n]  = size(cdat);
             cdati  = zeros(newns,n);
             disp(['Resampling section ' num2str(i) '  of  ' num2str(lf2c)]);
             for jj = 1:n
                 cdati(:,jj) = interp1(oldzz,cdat(:,jj),newzz);
             end
             dat{1,i} = cdati;
             zz{1,i} = newzz;
             ns(i) = newns;
             dzz(i) = newdz;
         end
     end
     % Check consistency of depth vector lengths and trim or pad with zeros
     if ~isempty(find(diff(ns)~=0)), %#ok<EFIND>
         oops = menu({'At least two sections have different', 'depth vectors!'}, ...
             'Zero-pad all to maximum depth','Trim all to minimum depth',...
             'Abort');
         if oops >= 3 || oops <= 0,
             disp('Operation Aborted!')
             return
         elseif oops == 1,
             newns = max(ns);
             disp('Padding with zeros')
         elseif oops == 2
             newns = min(ns);
         end
         for i=1:lf2c
             cdat    = dat{1,i};
             [m, n]  = size(cdat);
             cdati   = zeros(newns, n);
             if oops == 1,
                 cdati(1:m, 1:n) = cdat;
             elseif oops == 2
                 cdati(1:newns, 1:n) = cdat(1:newns, 1:n);
             end
             dat{1,i}    = cdati;
             ns(i) = newns;
             zs  = zz{1,i};
             zsi = zeros(1,newns);
             if oops == 1,
                 zsi(1:m) = zs;
             elseif oops == 2,
                 zsi(1:newns) = zs(1:newns);
             end
             zz{1,i} = zsi;
             disp(['Padding section ' num2str(i) '  of  ' num2str(lf2c)]);
         end
     end
 end
% Check consistency of trace spacings and resample if necessary
    if ~isempty(find(diff(dx)~=0)), %#ok<EFIND>
        disp('At least 2 sections have different trace spacings! ');
        disp('Resampling - new spacing will be the median of trace spacings')
        newdx = median(dx);
        for i=1:lf2c
           % Resample radargram first
            oldx    = scanax{1,i}; 
            minx    = oldx(1,1);
            maxx    = oldx(1,ntrs(i)); 
            newx    = minx:newdx:maxx; 
            newntr  = length(newx); 
            cdat    = dat{1,i};
            m   = size(cdat,1);
            cdati   = zeros(m,newntr);
            disp(['Re-spacing section = ' num2str(i) '  of  ' num2str(lf2c)]);
            for jj = 1:m
                cdati(jj,:) = interp1(oldx, cdat(jj,:), newx);
            end
            dat{1,i}    = cdati;
            scanax{1,i} = newx;
            ntrs(i)     = newntr;
            dx(i)       = newdx;
            % Resample trace x-data next
            xd  = xdata{1,i};
            xdi = interp1(oldx, xd, newx);
            xdata{1,i} = xdi;
        end
    end
% Check consistency of scan-axis lengths and pad with zeros
    if ~isempty(find(diff(ntrs)~=0)), %#ok<EFIND>
        oops = menu({'At least two sections have different', 'scan axis lengths!'}, ...
            'Zero-pad all to longest section','Trim all to shortest section (Recommended)',...
            'Abort');
        if oops >= 3 || oops <= 0,
            disp('Operation Aborted!')
            return
        elseif oops == 1,
            newntr = max(ntrs);
            disp('Padding with zeros')
        elseif oops == 2
            newntr = min(ntrs);
        end
        for i=1:lf2c
            cdat    = dat{1,i};
            [m, n]  = size(cdat);
            cdati   = zeros(m, newntr);
            if oops == 1,
                cdati(1:m, 1:n) = cdat;
            elseif oops == 2
                cdati(1:m, 1:newntr) = cdat(1:m, 1:newntr);
            end
            dat{1,i}    = cdati;
            ntrs(i) = newntr;
            xd  = xdata{1,i};
            xdi = zeros(newntr,1);
            if oops == 1,
                xdi(1:n) = xd;
                for j = n+1:newntr     
                    xdi(j) = xdi(j-1) + dx(i);
                end
            elseif oops == 2,
                xdi(1:newntr) = xd(1:newntr);
            end
            xdata{1,i} = xdi;
            disp(['Padding section = ' num2str(i) '  of  ' num2str(lf2c)]);
        end
    end
 % Check data for uniformity!
    if ~isempty(find(diff(tro)~=0)), %#ok<EFIND>
        erh = warndlg(['At least 2 sections have different antenna '...
            'separations! Proceed at your own risk!'], ...
            'MAKE3DDATA: WARNING');
        uiwait(erh);
    end
    % ALL is OK - create the data volume
    d3d = zeros(ns(1),ntrs(1),lf2c);
    yy = zeros(1,lf2c);
    % Sort sections in ascending-y order 
    for i=1:lf2c,
        yy(i) = ydata{1,i}(1);
    end
    [yy, iy] = sort(yy);
    % Create data volume in ascending-y order
    for i=1:lf2c,
        d3d(:,:,i) = dat{1,iy(i)};
    end
    d3d = permute(d3d,[3 2 1]);
    if isempty(find(dzz)), %#ok<EFIND>
        [x,y,z] = meshgrid(xdata{1}(:), yy, twt{1} );
        zlab = 'Time (ns)';
    end
    if isempty(find(dtt)), %#ok<EFIND>
        [x,y,z] = meshgrid(xdata{1}(:), yy, zz{1} );
        zlab = 'Depth (m)';
    end
% Export volume as a data structure
    D3D = struct('x',x, 'y',y, 'z',z, 'd', d3d, 'zlab', zlab);
% Wrap up and exit
    delete(findobj('tag','make3d_gui')); 
    return
end

function filebx(filebox,eventdata) %#ok<INUSD>
global files2concat
xfname = get(findobj('tag','filebox'),'string'); 
if exist(xfname)==2,  %#ok<EXIST>
    files2concat = get(findobj('tag','lsbox'),'string'); 
    nfiles = length(files2concat); 
    files2concat(nfiles+1) = cellstr(xfname);  
    set(findobj('tag','lsbox'),'string',files2concat); 
else
    erw = warndlg('No such file!','Warning'); 
    uiwait(erw); 
end
return

function brse(browsebtn,eventdata) %#ok<INUSD>
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

function rmfilefromlist(rmlsbutn,eventdata) %#ok<INUSD>
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

