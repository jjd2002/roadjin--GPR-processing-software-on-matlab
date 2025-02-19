function [v2d, zv, dzmig] = get2dvelocity(ns, ntr, xmax, dx)

%
% GET2DVELOCITY : Imports (from disk file) a synthetic structural model
%                 prepared by function BUILD2DMODEL and returns the
%                 corresponding velocity structure. The return velocity
%                 variable is a cell array containing:
%                 1. The non-dispersive term of the 2-D phase velocity
%                    model with values expressed in m/ns, computed after
%                    Bano (1996, Geophys. J. Int., 279 - 288).   
%                 2. The 2-D quality factor of the model structure. 
%                 3. Antenna central frequency (a scalar).
%               * The antenna frequency and the quality factor will be
%                 useful in incorporating the frequency dependence of the
%                 phase velocity (dispersion) in 2-D Split-step Fourier
%                 depth migration, (and future implementations of other 2-D
%                 spectral domain migration methods), also in the sense of   
%                 Bano (1996).
%
%         Usage : [v2d, zv, dzmig] = get2dvelocity(ns, ntr, xmax, dx)
%
%        Inputs : 
%            ns : The vertical dimension of the GPR section, equal to the
%                 number of samples per trace. 
%           ntr : The dimension of the scan axis, equal to the number of
%                 traces. 
%          xmax : Total length of the scan-axis in metres
%            dx : The trace spacing in metres.                          
%
%       Outputs : 
%           v2d : {1 x 3} Cell array containing: 
%                 1. v2d{1} : [ns x ntr] array, the non-dispersive term of
%                             the velocity model. 
%                 2. v2d{2} : [ns x ntr] array, the quality factor of the
%                             the 2-D model.  
%                 3. v2d{3} : Scalar, float, antenna central frequency.
%        zv(ns) : The z-coordinate axis of the velocity model.
%         dzmig : The depth stepping to be used for depth migration.
%
%      Requires : isinpoly.m, checkcomma.m
%
%      Author   : Andreas Tzanis
%                 Department of Geophysics, 
%                 University of Athens
%                 atzanis@geol.uoa.gr
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%

global ENVAR
% If the program is running alone, ENVAR was not initialized. Must define
% the endian!  
if isempty(ENVAR) || ~isfield(ENVAR,'endian'),
    [computer_system, maxsize, endian] = computer;
    endian = lower(endian);
    clear computer_system maxsize;
else
    endian = ENVAR.endian;
end
%%%%%   Ask what to do 
importmodel = questdlg('Import or Create Velocity Structure?',...
    'GET2DVELOCITY : REQUEST','Import','Create','Cancel','Create');
if strcmp(importmodel,'Cancel'),                   % Oops, abort!
    v2d = [];    zv = [];    dzmig = [];
    return
end

%%%%%   Velocity structure will be read from file
if strcmp(importmodel,'Import'),
    [vfname, vpname]= uigetfile('*.dat',...
        'GET2DVELOCITY : 2-D velocity model file name', ...
        ENVAR.currentworkdir);
    if isequal(vfname,0),
        v2d = [];    zv = [];    dzmig = [];
        return
    end
    fid    = fopen([vpname vfname],'r',endian);
    nsv    = fread(fid,1,'int');
    ntrv   = fread(fid,1,'int');
    antenna_fc = fread(fid,1,'float');
    dzmig  = fread(fid,1,'float');
    zv     = fread(fid,nsv,'float');
    vmod   = fread(fid,[nsv ntrv],'float'); 
    Q      = fread(fid,[nsv ntrv],'float');
    fclose(fid);
    v2d = { vmod Q antenna_fc };
    return
end

%%%%%   Velocity structure will be built from model file
nz   = ns;
m0   = 4*pi*1.0e-7;
e0   = 8.8592e-12;
bodies  = cell(1,4);            % initialize cell arrays to hold body data

% Get model file
[mfname, mpname]= uigetfile('*.dat','Give model file name', ...
    ENVAR.currentworkdir);
if mfname == 0, 
    v2d   = [];    zv = [];    dzmig = [];
    return, 
end
fid = fopen([mpname mfname],'r');
antenna_fc = fscanf(fid,' %f \n',1);      % get antenna frequency (in MHz)
dummy = fscanf(fid,' %f %f %f %f \n',4);  % read upper, lower model limits
xl = [dummy(1) dummy(2)];
zl = [dummy(3) dummy(4)];
% Load the model onto cell arrays
nbody = fscanf(fid,' %d ',1);             % read # bodies
for ib=1:nbody,                           % scan through bodies
    dummy = fgetl(fid);                   % read body tag
    bodies{ib,4} = double(dummy);
    dummy = fscanf(fid,' %d %f %f %f %f \n',5);
    nv    = dummy(1);                     % # vertices in this body
    bodies{ib,3} = dummy(2:5)';           % store EM properties
    x = zeros(nv,1);   z = zeros(nv,1);
    for iv = 1:nv,
        dummy = fscanf(fid,' %f  %f \n',2);
        x(iv)  = dummy(1);                % read vertices
        if x(iv) == xl(2),
            x(iv) = xmax;
        end
        z(iv)  = dummy(2);
    end
    bodies{ib,1} = x;
    bodies{ib,2} = z;
end
if xl(2) ~= xmax,
    xl(2) = xmax; %#ok<NASGU>
end
fclose(fid);

%%% Confirm antenna frequency to be used for this model
cldef  = cellstr(num2str(antenna_fc));
answer = inputdlg('Please CONFIRM ANTENNA FREQUENCY for this model', ...
    'GET2DVELOCITY : REQUEST',1 , cldef); 
if isempty(answer)
    v2d = [];    zv  = [];    dzmig = []; 
    return      
end
lb           = checkcomma(answer);     
antenna_fc = str2num(lb(1,:));
fc2          = 1.5*antenna_fc;

%%% Velocity grid size is fixed to ns x ntr
props = [];
for ib = 1:nbody
    props = [props; bodies{ib,3}]; %#ok<AGROW>
end
% test velocity
Vtest = 1.0./(sqrt(m0.*props(:,3).*e0.*props(:,1)).*cos((pi/4).*...
    (1.0-((2/pi).*atan(props(:,2))))));
Vtest = mean(Vtest);
% Determine possible grid size and spacings based on max properties and 
% antenna frequency
dzmig = Vtest*(1/(fc2*1e6))/8;
%%% Confirm depth-step size 
cldef  = cellstr(num2str(dzmig));
answer = inputdlg('Please CONFIRM STEP SIZE for depth migrations', ...
    'GET2DVELOCITY : REQUEST',1 , cldef); 
if isempty(answer)
    v2d = [];    zv  = [];    dzmig = []; 
    return      
end
lb     = checkcomma(answer);     
dzmig  = str2num(lb(1,:));

% Determine possible time sampling intervals for traveltime axis
dz  = (zl(2) - zl(1))/(nz-1);
zv  = 0:dz:(nz-1)*dz;
K   = zeros(nz,ntr);
Q   = zeros(nz,ntr);
M   = zeros(nz,ntr);

% - ------- Scan through the model to build the velocity grid ------------- 
%%% Set up x- and z- coordinate grids for "inpolygon" and "isinpoly"
[xtest, ztest] = meshgrid(0:dx:(ntr-1)*dx, 0:dz:(nz-1)*dz);
%%% The MATLAB function "inpolygon" may not exist for MATLAB releases
%%% earlier that 13! Check and if not, program will switch to Kirill
%%% Pankratov's "isinpoly".  
find_inpolygon = which('inpolygon');
if isempty(find_inpolygon),
    disp('GET2DVELOCITY > INPOLYGON not found! Working with ISINPOLY');
end
%%% Begin    
for ib=1:nbody,    
    x = bodies{ib,1};
    z = bodies{ib,2};
    props = bodies{ib,3};           % EM properties
    Kb = props(1);
    Qb = props(2);
    Mb = props(3);
    str = ['Scanning body ' num2str(ib) ' of ' num2str(nbody) '.'];
    msg = msgbox([str ' Please wait!'],...
          'Generating velocity profile','help');
%%% Fill with parameters of ib'th object. 
    if ~isempty(find_inpolygon),
        IN = inpolygon(xtest,ztest,x,z);
    elseif isempty(find_inpolygon),
        IN = isinpoly(xtest,ztest,x,z);
    end
    for i = 1:ns
        ii = find(IN(i,:));
        K(i,ii) = IN(i,ii)*Kb;
        Q(i,ii) = IN(i,ii)*Qb;
        M(i,ii) = IN(i,ii)*Mb;
    end;
    close(msg); 
end
% Model velocity profile (non-dispersive term)
vmod = 1.0./(sqrt(m0.*M*e0.*K).*cos((pi/4).*(1.0-((2/pi).*atan(Q)))))/1e9;
% Create output cell array
v2d = { vmod Q antenna_fc };

%%% Preview 
if exist('imagedisplay.m','file'),
    figure('Name','2-D Velocity Model', 'Menubar','none',...
        'NumberTitle','off', 'Tag','v2dfig');
    imagedisplay(vmod,0:dx:(ntr-1)*dx,zv,'Scanline (m)', 'Depth (m)');
end

%%%%%   Save velocity model
answer = questdlg('Save velocity model ?','GET2DVELOCITY : REQUEST','Yes');
if strcmp(answer,'Yes'),
    [vfname, vpname]= uiputfile('*.dat',...
        'GET2DVELOCITY : Save 2-D velocity model', ...
        ENVAR.currentworkdir);
    if isequal(vfname,0),
        return
    end
    fid = fopen([vpname vfname],'w',endian);
    fwrite(fid,nz,'int');
    fwrite(fid,ntr,'int');
    fwrite(fid,antenna_fc,'float');
    fwrite(fid,dzmig,'float');
    fwrite(fid,zv,'float');
    fwrite(fid,vmod,'float'); 
    fwrite(fid,Q,'float'); 
    fclose(fid);
end
delete(findobj('tag','v2dfig'));
return
