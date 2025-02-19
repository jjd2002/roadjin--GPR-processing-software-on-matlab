function [dmod,pos,tout,dx,dt] = run_FDTD_model(sig,ep,mu,x,z,t,Fc)
%
% RUN_FDTD_MODEL: Driver for the program TM_MODEL2D.M. 2-D, TM-mode, FDTD
% modeling for reflection GPR based on J.Inrving and R.Knight (2006,
% Computers and Geosciences, 32, 1247–1258). This program accepts matrices
% of the EM properties of a synthetic structural model (normally prepared
% by GETEMPROPERTIES.M), modifies them fr compliance with TM_MODEL2D and
% executes TM_MODEL2d.M.  
% ==> This program is part of the MATGPR suite.
%
% Usage  : [dmod,pos,tout,dx,dt] = run_FDTD_model(sig,ep,mu,x,z,t,Fc);
%
%  Input : sig  : Matrix of the model conductivity distribution with
%                 dimensions (length(z) x length(x)) 
%          ep   : Matrix of the model relative dielectric constant
%                 dimensions (length(z) x length(x))
%          mu   : Matrix of the model magnetic permeability with dimensions
%                 (length(z) x length(x))
%          x    : The x-axis of the synthetic model in metres
%          z    : The depth (z-axis) of the synthetic model in metres
%          t    : Initial estimate of the traveltime vector, to be used for
%                 determining the time step and number of iterations in the
%                 model radargram (in s) 
%          Fc   : The central frequency of the GPR antenna.
%
% Output : dmod : The synthetic (model) radargram
%          pos  : The vector of trace coordinates along the scan axis of
%                 the synthetic model 
%          tout : The 2-way traveltime of the synthetic model
%          dx   : Trace spacing of the synthetic model
%          dt   : The sampling rate of tt2w
%
%Requires: checkcomma.m, blackharrispulse.m, finddx.m, finddt.m,
%          gridinterp.m, padgrid.m  
%
% Author : Andreas Tzanis
%          Department of Geophysics, 
%          University of Athens
%          atzanis@geol.uoa.gr
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

% Irving and Knight's FDTD code wants the model matrices to be transposed
sig = sig';
ep  = ep';
mu  = mu';

% calculate minimum and maximum relative permittivity and permeability
% in the model (to be used in finddx.m and finddt.m) 
epmin = min(min(ep));
epmax = max(max(ep));
mumin = min(min(mu));
mumax = max(max(mu));

% create source pulse vector for use with finddx.m
srcpulse=blackharrispulse(Fc,t);

% use finddx.m to determine maximum possible spatial field discretization
% (in order to avoid numerical dispersion)
[dx,wlmin,fmax] = finddx(epmax,mumax,srcpulse,t,0.02);
disp(' ');
disp(['RUN_FDTD_MODEL > Maximum frequency contained in source pulse = ',...
    num2str(fmax/1e6),' MHz']);
disp(['RUN_FDTD_MODEL > Minimum wavelength in simulation grid = ',...
    num2str(wlmin),' m']);
disp(['RUN_FDTD_MODEL > Maximum possible electric/magnetic field ', ... 
    'discretization (dx,dz) = ',num2str(dx),' m']);
disp(['RUN_FDTD_MODEL > Maximum possible electrical property ', ... 
    'discretization (dx/2,dz/2) = ',num2str(dx/2),' m']);
disp(' ');
% Choose model discretization
clb(1) = cellstr(['Give trace spacing (dx). Maximum Possible is ',num2str(dx),' m']);
clb(2) = cellstr(['Give depth spacing (dz). Maximum Possible is ',num2str(dx/2),' m']);
clb(3) = cellstr(['Give source start location, 0<=x0<=',num2str(max(x)),' m']);
clb(4) = cellstr(['Give source Spacing. Min = dx, Max = ',num2str(max(x)/2),' m']);
cldef(1) = cellstr(num2str( fix(1000*dx)/1000 ));   
cldef(2) = cellstr(num2str( fix(1000*dx/2)/1000 ));
cldef(3) = cellstr('0');
cldef(4) = cellstr(num2str( fix(1000*dx*2)/1000 ));
answer = inputdlg(clb,'Give Final Model Discretization',1,cldef );
if isempty(answer), 
    dmod=[]; pos=[]; tout=[]; dx=[]; dt=[];
    return
end
lb    = checkcomma(answer);
dx    = str2num(lb(1,:));
dz    = str2num(lb(2,:));
x0    = str2num(lb(3,:));
dxsrc = str2num(lb(4,:));
disp(['RUN_FDTD_MODEL > Using dx = ',num2str(dx),' m,  dz = ',...
    num2str(dz),' m   x0 = ',num2str(x0),' m  dxsrc = ',num2str(dxsrc)]);

% find the maximum possible time step using this dx and dz
% (in order to avoid numerical instability)
tmax     = max(t);
dtmax    = finddt(epmin,mumin,dx,dz);
disp(['RUN_FDTD_MODEL > Maximum possible time step with this ', ...
    'discretization = ', num2str(dtmax/1e-9),' ns']);
disp(' ');
% Choose time step
clb      = cell(1,2);
clb(1)   = cellstr(['Give sampling interval dt. Maximum Possible is ',...
                     num2str(dtmax/1e-9),' ns']);
clb(2)   = cellstr('Give total 2-way traveltime');
cldef(1) = cellstr(num2str( fix(1000*(dtmax/1e-9))/1000 ));   
cldef(2) = cellstr(num2str( tmax/1e-9 ) );   
answer   = inputdlg(clb,'Give Time Step',1,cldef );
if isempty(answer), 
    dmod=[]; pos=[]; tout=[]; dx=[]; dt=[];
    return
end
lb       = checkcomma(answer);
dt       = str2num(lb(1,:))*1.0e-9;
tmax     = str2num(lb(2,:))*1.0e-9;
disp(['RUN_FDTD_MODEL > Using dt = ',num2str(dt/1e-9),' ns,  tmax = ',...
    num2str(tmax/1e-9),' ns']);
disp(' ');
% create final time vector (s) and corresponding source pulse
% (using the proper values of dt and tmax this time)
t        = 0:dt:tmax;
srcpulse = blackharrispulse(Fc,t);    

% interpolate electrical property grids to proper spatial discretization
% NOTE:  we MUST use dx/2 here because we're dealing with electrical
% property matrices 
disp('RUN_FDTD_MODEL > Interpolating electrical property matrices...');
disp(' ');
x2 = min(x):dx:max(x);  
% Make sure that dimension of x2 is odd 
if mod(max(size(x2)),2)~=1,
    nx = max(size(x2)) - 1;
    dx = (max(x) - min(x))/(nx-1);
    x2 = min(x):dx:max(x); 
end
z2 = min(z):dz:max(z);
% Make sure that dimension of z2 is odd 
if mod(max(size(z2)),2)~=1,
    nz = max(size(z2)) - 1;
    dz = (max(z) - min(z))/(nz-1);
    z2 = min(z):dz:max(z); 
end
ep2  = gridinterp(ep,x,z,x2,z2,'nearest');
mu2  = gridinterp(mu,x,z,x2,z2,'nearest');
sig2 = gridinterp(sig,x,z,x2,z2,'nearest');

% pad electrical property matrices for PML absorbing boundaries
npml = 10;  % number of PML boundary layers
ep3          = padgrid(ep2,x2,z2,2*npml+1);
mu3          = padgrid(mu2,x2,z2,2*npml+1);
[sig3,x3,z3] = padgrid(sig2,x2,z2,2*npml+1);

% create source and receiver location matrices
% (rows are [x location (m), z location (m)])
srcx = (x0:dxsrc:max(x))';
srcz = 0*ones(size(srcx));
recx = srcx;                                   % Assume Tx and Rx coincide!
recz = srcz;
srcloc = [srcx srcz];
recloc = [recx recz];

% clear unnecessary matrices taking up memory
 clear x x2 z z2 ep ep2 mu mu2 sig sig2 

% Output and Plotting options
clb(1) = cellstr(['Write a sample to the output matrix every # iterations']);
clb(2) = cellstr(['View in real time - 0=NO, 1=YES']);
cldef(1) = cellstr('1');   
cldef(2) = cellstr('1');   
answer = inputdlg(clb,'Output and Plotting Options',1,cldef );
if isempty(answer), 
    dmod=[]; pos=[]; tout=[]; dx=[]; dt=[];
    return
end
lb           = checkcomma(answer);
outstep      = str2num(lb(1,:));
plotopt(1)   = str2num(lb(2,:));

% Set Plotting defaults
plotopt(2)    = 50;
plotopt(3)    = 0.05;
movieopt.ans  = 0;                                    % Default is no movie
% or change plotting options
if plotopt(1) == 1,
    clb(1) = cellstr(['Plot every # iterations']);
    clb(2) = cellstr(['Colour Saturation Threshold']);
    clb(3) = cellstr(['Make Movie, 0=NO, 1=YES ']);
    cldef(1) = cellstr(num2str(plotopt(2)));
    cldef(2) = cellstr(num2str(plotopt(3)));
    cldef(3) = cellstr(num2str(movieopt.ans));
    answer = inputdlg(clb,'Plotting Options',1,cldef );
    if isempty(answer),
        dmod=[]; pos=[]; tout=[]; dx=[]; dt=[];
        return
    end
    lb           = checkcomma(answer);
    plotopt(2)   = str2num(lb(1,:));
    plotopt(3)   = str2num(lb(2,:));
    movieopt.ans = str2num(lb(3,:));

    if movieopt.ans == 1,
        clb    = cell(4,1);
        clb(1) = cellstr('AVI file name');
        clb(2) = cellstr('AVI Compression - 1 = "Cinepak", 0=None');
        clb(3) = cellstr('Frames per second');
        clb(4) = cellstr('Movie Quality');
        cldef(1) = cellstr('work\Ey_wavefield.avi');
        cldef(2) = cellstr('1');
        cldef(3) = cellstr('2');
        cldef(4) = cellstr('80');
        answer = inputdlg(clb,'Movie parameters',1,cldef );
        if isempty(answer),
            dmod=[]; pos=[]; tout=[]; dx=[]; dt=[];
            return
        end
        lb             = checkcomma(answer);
        movieopt.fname = lb(1,:);
        cpropt = str2num(lb(2,:));
        if cpropt,
            movieopt.cpr = 'Cinepak';
        else
            movieopt.cpr = 'None';
        end
        movieopt.fps   = str2num(lb(3,:));
        movieopt.q     = str2num(lb(4,:));
    end
end

% run the simulation
tic;
[gather,tout,srcx,srcz,recx,recz] = ...
    TM_model2d(ep3,mu3,sig3,x3,z3,srcloc,recloc,srcpulse,t,npml,outstep,...
    plotopt, movieopt);
disp(' ');
disp(['RUN_FDTD_MODEL > Total running time = ',num2str(toc/3600),' hours']);

% extract common offset reflection GPR data from multi-offset data cube 
dmod = zeros(max(size(tout)), size(srcx,1));
for i=1:length(srcx);
    dmod(:,i) = gather(:,i,i);
end
pos = (srcx+recx)/2;

return