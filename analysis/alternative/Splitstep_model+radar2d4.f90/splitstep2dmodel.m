function [dm,scanaxis,tt2w,nx,dx,nt,dt,antenna_freq] = splitstep2dmodel()
%
% SPLITSTEP2DMODEL: GPR forward modelling using the method of Bitri and
% Grandjean (1998,Geophysical Prospecting, 46, 287-301). This program
% imports (from disk file) a synthetic structural model prepared by
% function build2dmodel.m and generates matrices of the relative dielectric
% constant, magetic permeability and quality factor structures. These are
% passed to a Fortran-90 MEX-file that computes the corresponding velocity
% structure after Bano (1996, Geophys. J. Int., 279 - 288) and hence the
% forward model. 
% ==> This program is part of the MATGPR suite.
%
% Usage  : [dm,scanaxis,tt2w,nx,dx,nt,dt,antenna_freq] = splitstep2dmodel;
%
% Output : dm           : The synthetic GPR section 
%          tt2w         : The 2-way traveltime of the synthetic model
%          nt           : The dimension of tt2w  
%          dt           : The sampling rate of tt2w
%          scanaxis     : The vector of traces coordinates along the scan
%                         axis of the synthetic model
%          nx           : The dimension of scanaxis
%          dx           : Trace spacing of the synthetic model
%          antenna_freq : The central frequency of the GPR antenna.
%
%Requires: checkcomma.m, isinpoly.m (for MATLAB R12 and earlier)
%
% Author : Andreas Tzanis
%          Department of Geophysics, 
%          University of Athens
%          atzanis@geol.uoa.gr
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

m0   = 4*pi*1.0e-7;
e0   = 8.8592e-12;
bodies  = cell(1,4);            % initialize cell arrays to hold body data
% Get model file
[mfname, mpname]= uigetfile('*.dat','Give model file name');
if mfname == 0, 
    dm = [];  scanaxis = [];  tt2w = [];
    nt = [];  dt = [];  nx = [];  dx = [];  antenna_freq = [];
    return, 
end
fid = fopen([mpname mfname],'r');
antenna_freq = fscanf(fid,' %f \n',1);     % get antenna frequency (in MHz)
fc2          = 1.5*antenna_freq;
dummy = fscanf(fid,' %f %f %f %f \n',4);   % read upper, lower model limits
xl = [dummy(1) dummy(2)];
zl = [dummy(3) dummy(4)];
% Load the model onto cell arrays
nbody = fscanf(fid,' %d ',1);                     % read # bodies
for ib=1:nbody,                                   % scan through bodies
    x = [];
    z = [];
%    dummy = fscanf(fid,' %s \n',1)               % read body tag
    dummy = fgetl(fid);
    bodies{ib,4} = double(dummy);
    dummy = fscanf(fid,' %d %f %f %f %f \n',5);
    nv    = dummy(1);                             % # vertices in this body
    bodies{ib,3} = dummy(2:5)';                   % EM properties
    for iv = 1:nv,
        dummy = fscanf(fid,' %f  %f \n',2);
        x(iv)  = dummy(1);                        % read vertices
        z(iv)  = dummy(2);
    end
    bodies{ib,1} = x';
    bodies{ib,2} = z';
end
fclose(fid);

% Compute an appropriate grid size based on max, average or minimum 
% model properties
props = [];
for ib = 1:nbody
    props = [props; bodies{ib,3}];
end
% test velocity
vv = 1.0./(sqrt(m0.*props(:,3).*e0.*props(:,1)).*cos((pi/4).*...
    (1.0-((2/pi).*atan(props(:,2))))));
Vtest = [max(vv) mean(vv) min(vv)];
% Determine possible grid sizes and spacings 
dz = Vtest*(1/(fc2*1e6))/8;
nz = ceil( (zl(2) - zl(1))./dz );
dx = 2.0*dz;
nx = ceil( (xl(2) - xl(1))./dx );      
% Determine possible sampling rates for the traveltime axis
dt=1/(1.5*(fc2*1.0e6)*4);
totalt = (zl(2)-zl(1))*2.0./Vtest;       % 2-way traveltime
nt = ceil(1.1*totalt/dt);                % Size of time / frequency axis
dt = dt*1e9;                             % Scale dt for passing to RADAR2D
% Choose model size
clb(1) = cellstr('Suggested velocity grid size: ');
clb(2) = cellstr('     LARGE      MEDIUM       SMALL   ');
clb(3) = cellstr(['      ' num2str(nz(3)) 'x' num2str(nx(3)) ... 
                  '         ' num2str(nz(2)) 'x' num2str(nx(2)) ...
                  '         ' num2str(nz(1)) 'x' num2str(nx(1)) ]);
clb(4) = cellstr('Suggested model radargram size: ');
clb(5) = cellstr('     LARGE      MEDIUM       SMALL  ');
clb(6) = cellstr(['      ' num2str(nt(3)) 'x' num2str(nx(3)) ... 
                  '         ' num2str(nt(2)) 'x' num2str(nx(2)) ...
                  '         ' num2str(nt(1)) 'x' num2str(nx(1)) ]);
np = menu(clb,'Small','Medium','Large','CUSTOM','Cancel');    
switch np
case 3
    np = 3;
    dz = dz(np);  nz = nz(np);  dx = dx(np);  nx = nx(np);  nt = nt(np);
case 2
    np = 2;
    dz = dz(np);  nz = nz(np);  dx = dx(np);  nx = nx(np);  nt = nt(np);
case 1
    np = 1;
    dz = dz(np);  nz = nz(np);  dx = dx(np);  nx = nx(np);  nt = nt(np);
case 4                                           % User's own choise
    clb = cell(4,1);
    clb(1) = cellstr('Give model X-dimension (nx)');
    clb(2) = cellstr('Give model Z-dimension (nz)');
    clb(3) = cellstr('Give traveltime dimension (nt)');
    clb(4) = cellstr('Give sampling rate (dt)');
    cldef(1) = cellstr(num2str(nx(2)));          % defaults are mean values
    cldef(2) = cellstr(num2str(nz(2)));                            
    cldef(3) = cellstr(num2str(nt(2)));                            
    cldef(4) = cellstr(num2str(dt));                            
    answer = inputdlg(clb,'Give Model Grid Sizes',1,cldef ); 
    if isempty(answer),  
        dm = [];  scanaxis = [];  tt2w = [];
        return, 
    end;
    lb = checkcomma(answer); 
    nx  = str2num(lb(1,:));      
    dx  = (xl(2) - xl(1))/(nx-1);
    nz  = str2num(lb(2,:));
    dz  = (zl(2) - zl(1))/(nz-1);
    nt  = str2num(lb(3,:));
    dt  = str2num(lb(4,:));

case 5                                           %case 'Cancel'
    dm = [];  scanaxis = [];  tt2w = [];
    nt = [];  dt = [];  nx = [];  dx = [];  antenna_freq = [];
    return
end

tic                                              % start clock
% Initialize grids of model parameters 
K = zeros(nz,nx);
Q = zeros(nz,nx);
M = zeros(nz,nx);

%%% ------- Scan through the model to build the velocity grid ------------- 
%%% Set up x- and z- coordinate grids for "inpolygon" and "isinpoly"
[xtest, ztest] = meshgrid(xl(1):dx:(nx-1)*dx, zl(1):dz:(nz-1)*dz);
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
    props = bodies{ib,3};                   % EM properties
    Kb = props(1);
    Qb = props(2);
    Mb = props(3);
    str = ['Scanning object ' num2str(ib) ' of ' num2str(nbody) '.'];
    msg = msgbox([str ' Please wait!'],...
          'Generating velocity profile','help');
%%% Fill with parameters of ib'th object. 
    if ~isempty(find_inpolygon),
        IN = inpolygon(xtest,ztest,x,z);
    elseif isempty(find_inpolygon),
        IN = isinpoly(xtest,ztest,x,z);
    end
    for i = 1:nz
        ii = find(IN(i,:));
        K(i,ii) = IN(i,ii)*Kb;
        Q(i,ii) = IN(i,ii)*Qb;
        M(i,ii) = IN(i,ii)*Mb;
    end;
    close(msg); 
end
time = toc;
disp(['> K, Mu and Q grids computed in ' num2str(time) ' seconds'])

% Prepare distance and traveltime vectors
scanaxis   = xl(1):dx:(nx-1)*dx;
tt2w       = 0:dt:(nt-1)*dt;

%%%% Compute model 
disp('> Beginning model computations')
tic
% Create model files for RADAR2D4
efile   = tempname;
fid=fopen(efile,'w','l');
for i=1:nx; 
    fwrite(fid,K(:,i),'float'); 
end
fclose(fid);
qfile   = tempname;
fid=fopen(qfile,'w','l');
for i=1:nx; 
    fwrite(fid,Q(:,i),'float'); 
end
fclose(fid);
mgfile  = tempname;
fid=fopen(mgfile,'w','l');
for i=1:nx; 
    fwrite(fid,M(:,i),'float'); 
end
fclose(fid);

% Run Radar2d4 
resfile = tempname;                  % results file
argstr = [efile ' ' mgfile ' ' qfile ' ' resfile ' ' num2str(nx) ' ' ...
    num2str(nt) ' ' num2str(dx,4) ' ' num2str(nz) ' ' ... 
    num2str(dz,4) ' ' num2str(dt,4) ' ' num2str(antenna_freq) ' 1'];
disp(['Radar2d4 arguments : ' argstr]);
radar2dpath = which('radar2d4.exe');
[busted, result] = dos([radar2dpath ' ' argstr],'-echo');
if busted,
    disp('SPLITSTEP2DMODEL > Ungraceful termination of RADAR2D4.EXE')
    disp('SPLITSTEP2DMODEL > due to system error.')
    dm = [];
    delete(efile);
    delete(qfile);
    delete(mgfile);   
    return
end
% Succesful termination of radar2d4.exe   
dm = load(resfile);                                % load results
dm = dm';

% Clear temporary work files
delete(efile);
delete(qfile);
delete(mgfile);
delete(resfile);
% wrap up
time  = toc;
disp(['> Model computed in ' num2str(time) ' seconds'])
return
