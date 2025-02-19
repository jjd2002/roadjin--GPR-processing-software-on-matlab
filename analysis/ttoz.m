function [zv,nz,dz,dm] = ttoz(traveltime,d,dt,vofh,firstt,firstz)
%
%   TTOZ  : Remap or resample a GPR section from traveltime vs. distance to
%           depth vs. distance, assuming a uniform or layered velocity
%           structure. The program:
%           1. Computes a time vs depth function t(z) from the velocity vs.
%              depth function v(t). 
%           2. Inverts t(z) to z(t).
%           3. Remaps z(t) to depth z and, based on this mapping,
%           4. Remaps E(t) --> E(z), where E stands for the amplitude
%              of a GPR trace.
%
%  Usage  : [zv,nz,dz,dm] = ttoz(traveltime,d,dt,vofh,firstt,firstz)
%
%  Inputs : 
%       d : [ns x ntr] matrix of the common-offset GPR section.
%      dt : Time sampling rate in nanosecs. 
%    vofh : [nlayer x 2 ] velocity model to be used for migration.
%           Two-column vector of "nlay" velocity - thickness pairs like:
%           [ 0.1    1 ;       ... 1st layer
%             0.08   2 ;       ... 2nd layer
%             ...
%             0.18   0 ]       ... n'th layer==basal halfspace
%           A uniform halfspace is given as a single layer stucture with
%           zero thickness. Velocity values in m/ns and thichnesses in m.
%  firstt : first time sample from which to begin migrating.
%  firstz : first z sample from which to begin migrating.
%
% Outputs :  
%      dm : [nz x ntr] matrix of the depth-migrated GPR section.
%      nz : The number of depth samples in "dm", computed internally. In
%           the current version of the program, nz = ns.
%      dz : The depth sampling interval in m.
%      zv : [1 x nz ] vector of the depth coordinates of "dm".
%      tz : times mapped onto depths zv
% 
%Requires : yxtoxy.m
%
%  Author : Andreas Tzanis
%           Department of Geophysics, 
%           University of Athens
%           atzanis@geol.uoa.gr
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

[ns,ntr] = size(d);
%%%%%   Velocity structure
[nlay, vhpairs] = size(vofh);
if vhpairs ~=2,
    erh = errordlg('Error in the structure of velocity model',...
        'TTOZ : ERROR');
    uiwait(erh)
    dm = [];   tz=[];   zv=[];   dz=[];
    return
end
layer_velocity  = vofh(:,1);
layer_thickness = vofh(:,2);
%%%%%   Check imported velocities
for i=1:nlay
    if layer_velocity(i) <= 0 | layer_velocity(i) > 0.2998,
        errordlg('Impossible velocity value found! Please try again!',...
            'TTOZ : ERROR')
        dm = [];   tz=[];   zv=[];   dz=[];
        return
    end
end
%%%%  Begin by converting given v(z) to t(z)
%%%%  Ensure an upper limit for t(z) 
if nlay == 1 & layer_thickness(nlay) == 0, 
    layer_thickness(nlay) = layer_velocity(nlay)*traveltime(ns);
else
    layer_thickness(nlay) = 2*sum(layer_thickness);
end
dz   =  (sum(layer_thickness)-firstz)/(ns - 1);
z    =  firstz:dz:sum(layer_thickness); 
nz   =  length(z);
tz(1) = 2.0*firstz/layer_velocity(1);
for iz=1:length(layer_thickness); 
    tz(iz+1) = tz(iz) + 2.0*layer_thickness(iz)/layer_velocity(iz);
end
tofz  = interp1([0; cumsum(layer_thickness)],tz',z);
%%%%%   Now compute z(t) from t(z)
vfz   = layer_velocity(1);
if nlay == 1,
    vlz   = vfz;
else
    vlz   = layer_velocity(length(layer_velocity)-1);
end
lt    = firstt+(ns-1)*dt;
lz    = firstz+(nz-1)*dz;
zt    = yxtoxy(nz,dz,firstz,tofz,ns,dt,firstt,0.0,0.0);
ii    = find(traveltime < tofz(1));
if ~isempty(ii),
    zt(ii) = 0.5*traveltime(ii)*vfz;
end
ii    = find(traveltime >= tofz(nz));
if ~isempty(ii),
    zt(ii) = lz + 0.5*(traveltime(ii) - tofz(nz))*vlz;
end
%%%%  Compute depth vector and then interpolate to t(z)
dz    = zt(ns)/(ns-1);
zv    = [firstz:dz:(ns-1)*dz]';
tz    = interp1(zt', traveltime', zv);

%%%%  Output array has same size as input 
%%%%  Note that there's a possibility of aliasing!
dm    = d;
nz    = ns;
dz    = zt(ns)/(nz-1);
%%%%  If dz is chosen on the basis of the minimum velocity, 
%%%%  aliasing may be suppressed
% dz  = min(layer_velocity)*dt/2.0;     
% nz = ceil(zt(ns)/dz) + 1;     
% dm = zeros(nz,ntr);

%%%%%%  Now loop over traces to map time to depth by resampling 
h = waitbar(0,'Resampling T -> Z');
for i=1:ntr
    tr = d(:,i);
    tr1 = interp1(traveltime,tr,tz,'cubic');
    dm(:,i) = tr1;
    waitbar(i/ntr,h);
end
close(h);
        
return
