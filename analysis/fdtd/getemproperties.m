function [S,K,M,scanaxis,zv,tt2w,Fc,nx,dx,nt,dt,nz,dz] = getemproperties()
%
% GETEMPROPERTIES : Imports (from disk file) a synthetic structural model
%                   prepared by function BUILD2DMODEL.M and discretizes it,
%                   constructing matrices of its EM properties. The
%                   dimensions of the matrices (discretized model) are
%                   determined automatically, with respect to 
%                   central frequency Fc of the probing GPR antenna.
%
%   Usage : [S,K,M,scanaxis,zv,tt2w,Fc,nx,dx,nt,dt,nz,dz] =
%                      getemproperties()
%
%  Outputs : 
%        S : Matrix of the model conductivity distribution with
%            dimensions (nz x nx)  
%        K : Matrix of the model relative dielectric constant with
%            dimensions (nz x nx)
%        M : Matrix of the model relative magnetic permeability with
%            dimensions (nz x nx)
% scanaxis : The x-axis of the synthetic model in metres
%       zv : The depth (z-axis) of the synthetic model in metres
%     tt2w : An initial estimate of the 2-way traveltime vector, to be used
%            for determining the time step and number of iterations for  
%            modelling (in seconds) 
%       Fc : The central frequency of the GPR antenna.
%       nx : The dimension of scanaxis
%       dx : Initial estimate of spacing (discretization step)along the
%            scanaxis in metres 
%       nt : The dimension of tt2w  
%       dt : Initial estimate of the sampling rate in tt2w (in seconds)
%       nz : The dimension of zv
%       dz : Initial estimate of vertical spacing (discretization step)
%            along zv in metres 
%
%      Requires : isinpoly.m, checkcomma.m
%
%      Author   : Andreas Tzanis
%                 Department of Geophysics, 
%                 University of Athens
%                 atzanis@geol.uoa.gr
%
% Copyright (C) 2008, Andreas Tzanis. All rights reserved.
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
m0   = 4*pi*1.0e-7;
e0   = 8.8592e-12;
bodies  = cell(1,4);            % initialize cell arrays to hold body data

% Get model file
[mfname, mpname]= uigetfile('*.dat','Give model file name', ...
    ENVAR.currentworkdir);
if mfname == 0, 
    K = [];  M = []; S = []; scanaxis = [];  tt2w = [];  zv = [];
    nt = [];  dt = [];  nx = [];  dx = [];  Fc = [];
    return, 
end
fid = fopen([mpname mfname],'r');
Fc = fscanf(fid,' %f \n',1);              % get antenna frequency (in MHz)
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
    bodies{ib,3} = dummy(2:5)';           % EM properties
    x = zeros(nv,1);   z = zeros(nv,1);
    for iv = 1:nv,
        dummy = fscanf(fid,' %f  %f \n',2);
        x(iv)  = dummy(1);                % read vertices
        z(iv)  = dummy(2);
    end
    bodies{ib,1} = x';
    bodies{ib,2} = z';
end
fclose(fid);

%%% Confirm antenna frequency to be used for this model
cldef  = cellstr(num2str(Fc));
answer = inputdlg('Please CONFIRM ANTENNA FREQUENCY for this model', ...
    'GETEMPROPERTIES',1 , cldef); 
if isempty(answer)
    K = [];  M = []; S = []; scanaxis = [];  tt2w = [];  zv = [];
    nt = [];  dt = [];  nx = [];  dx = [];  Fc = [];
    return      
end
lb   = checkcomma(answer);     
Fc   = str2num(lb(1,:));
fc2  = 1.5*Fc;

props = zeros(nbody, 4);
for ib = 1:nbody
    props(ib,:) = bodies{ib,3};
end
% Determine initial grid sizes and spacings 
% Choose good discretization to begin with
Vtest = 1.0./(sqrt(m0.*props(:,3).*e0.*props(:,1)).*cos((pi/4).*...
    (1.0-((2/pi).*atan(props(:,2))))));
dz = min(Vtest)*(1/(fc2*1e6))/8;
nz = ceil( (zl(2) - zl(1))./dz );
dx = 2.0*dz;
nx = ceil( (xl(2) - xl(1))./dx );      
% Determine initial sampling rate for the traveltime axis
dt=1/(1.5*(fc2*1.0e6)*4);
totalt = (zl(2)-zl(1))*2.0./median(Vtest);               % 2-way traveltime
nt = ceil(1.1*totalt/dt);                   % Size of time / frequency axis

% Make sure that the dimensions of the model grid are even
if mod(nx,2)==1,
    nx = nx-1;
    dx = (xl(2) - xl(1))/(nx-1);
end
if mod(nz,2)==1,
    nz = nz-1;
    dz = (zl(2) - zl(1))/(nz-1);
end

% Prepare distance, traveltime and depth vectors
Fc       = Fc*1.0e6;
scanaxis = xl(1):dx:(nx-1)*dx;
tt2w     = 0:dt:(nt-1)*dt;
zv       = 0:dz:(nz-1)*dz;
% Initialize grids of model parameters 
K   = zeros(nz,nx);
S   = zeros(nz,nx);
M   = zeros(nz,nx);
Q   = zeros(nz,nx);

% - ------- Scan through the model to build the velocity grid ------------- 
%%% Set up x- and z- coordinate grids for "inpolygon" and "isinpoly"
[xtest, ztest] = meshgrid(scanaxis, zv);
%%% The MATLAB function "inpolygon" may not exist for MATLAB releases
%%% earlier that 13! Check and if not, program will switch to Kirill
%%% Pankratov's "isinpoly".  
find_inpolygon = which('inpolygon');
if isempty(find_inpolygon),
    disp('GETEMPROPERTIES > INPOLYGON not found! Working with ISINPOLY');
end
%%% Begin    
for ib=1:nbody,    
    x = bodies{ib,1};
    z = bodies{ib,2};
    props = bodies{ib,3};           % EM properties
    Kb = props(1);
    Qb = props(2);
    Mb = props(3);
    Sb = props(4);
    str = ['Scanning body ' num2str(ib) ' of ' num2str(nbody) '.'];
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
        S(i,ii) = IN(i,ii)*(1/Sb);
        M(i,ii) = IN(i,ii)*Mb;
        Q(i,ii) = IN(i,ii)*Qb;
    end;
    close(msg); 
end

% Draw the model
preview_properties(S,K,M,Q,scanaxis,zv,Fc)
return

function preview_properties(S,K,M,Q,scanaxis,zv,Fc)

% Plots EM properties and velocity structure of a 2-D model

m0   = 4*pi*1.0e-7;
e0   = 8.8592e-12;
if isempty(Q),
    w     = 2*pi*Fc;
    sigmastar = S + sqrt(-1)*w*e0*K;
    estar     = -sqrt(-1)*sigmastar/w;
    Q = 1./(S -w*imag(estar))./(w*real(estar));
end
% Model velocity profile (non-dispersive term)
vmod = 1.0./(sqrt(m0.*M*e0.*K).*cos((pi/4).*(1.0-((2/pi).*atan(Q)))))/1e9;

% plot electrical property grids - use adaptive figure sizing and
% posistioning 
scrsz  = get(0,'screensize');
fpos = round([57*scrsz(3)/1680 525*scrsz(4)/1050 ...
    560*scrsz(3)/1680 420*scrsz(4)/1050]);
figure('Name','2-D Earth Model', 'Menubar','none',...
    'NumberTitle','off', 'Tag','emsfig', 'position', fpos);
figuretools;
imagecolors;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(221);
imagesc(scanaxis,zv,K);
xlabel('x (m)');   ylabel('Depth (m)');
title('\epsilon_r_e_l_a_t_i_v_e')
subplot(222);
imagesc(scanaxis,zv,log10(S));
xlabel('x (m)');   ylabel('Depth (m)');
title('log_1_0(\sigma) in S/m')
subplot(223);
imagesc(scanaxis,zv,M);
xlabel('x (m)');   ylabel('Depth (m)');
title('\mu_r_e_l_a_t_i_v_e')
subplot(224);
imagesc(scanaxis,zv,vmod); 
xlabel('x (m)');   ylabel('Depth (m)');
title('Velocity (nondispersive term in m/ns)');
return
