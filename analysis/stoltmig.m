 function dmig = stoltmig(d, dt, dx, vofh, interpm)
%
% STOLTMIG : Compact implementation of Stolt's F-K migration for uniform or
% layered velocity structures. The program uses Stolt stretching to
% accomodate non-uniform (layered) velocity structures, while the actual
% migration code was built with tips from Jon Claerbout's implementation of
% the Stolt algorithm (see Imaging the Earth's Interior, 1996, pp 55-57).
% However, in contrast to his sequential RATFOR programming, this
% implementation is fully vectorized and fast enough, soas to compete with
% compiled Fortran or C Stolt migration routines.
%
% Usage    : dmig = stoltmig(d, dt, dx, vofh, interpm)
%
% Inputs   : d         is the zero- or common-offset GPR section 
%          : dt        is the sampling rate 
%          : dx        is the trace spacing (spatial sampling rate) 
%          : vofh(n,2) is the 1-D velocity model of "nlay" velocity -
%                      thickness pairs, for example:
%                      [ 0.1    1 ;       ... 1st layer
%                        0.08   2 ;       ... 2nd layer
%                        ...
%                        0.18   0 ]       ... n'th layer==basal halfspace
%                      A uniform halfspace is given as a single layer
%                      stucture with zero thickness. Velocity values are
%                      given in m/ns and thichnesses in m.
%          : interpm   defines the interpolation method to be used. Any
%                      of the methods implemented in MATLAB's "interp1"
%                      routine can be used, i.e. 'nearest', 'linear',
%                      'spline', 'pchip', 'cubic' and 'v5cubic'. Default
%                      is the bandlimited linear interpolation (very fast
%                      and stable). Best choice is 'cubic' (slower).
%                      Least recommended are 'spline' and 'nearest'. 
%
% Output   : dmig      is the migrated section
%
% Requires : yxtoxy.m
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

% Check input arguments
if nargin < 5, 
    interpm = 'linear';       % default interpolation method 
end
if nargin < 4,
    erh = errordlg('Essential information missing! Please check input', ...
        'STOLTMIG : ERROR');
    uiwait(erh)
    dmig = []; 
    return
end

[ns,ntr]=size(d);
fn = 1/(2*dt);                     % Nyquist frequency
Range  =  (ns-1)*dt;               % total two-way traveltime
%%%%%   Velocity structure
[nlay, vhpairs] = size(vofh);
if vhpairs ~=2,
    erh = errordlg('Error in the structure of velocity model',...
        'STOLTMIG : ERROR');
    uiwait(erh)
    dmig = [];
    return
end
layer_velocity  = vofh(:,1)';
layer_thickness = vofh(:,2)';
for i=1:nlay
    if layer_velocity(i) <= 0 | layer_velocity(i) > 0.2998,
        errordlg('Impossible velocity value found! Please try again!',...
            'STOLTMIG : ERROR')
        dmig = []; 
        return
    end
end

%%%%%  Begin 
tic                                      % start clock
if nlay == 1,                            % Case of uniform halfspace
    vmig   = layer_velocity(1);
    dz = (vmig*Range/2)/(ns-1);          % vmig*Range/2 is max penetration!
    delt     = dt;
%%%%%  Pad and FFT
    disp('STOLTMIG > Transforming TX -> FKx');
    zpad = zeros(ns,ntr);
    fk   = fft2([d ; zpad]);  
    clear zpad;
    fk   = fftshift(fk);
    [nss,ntr] = size(fk);

elseif nlay > 1,                          % Case of layered halfspace
    firstt = 0.0;
    firstz = 0.0;
%%%%%   Compute RMS velocity from given V(z) function 
    zmax = max(vofh(:,1))*Range/2;        % max possible penetration
    delz = (zmax-firstz)/(ns-1);
    z    = firstz : delz : zmax;
    nz   = length(z);
    tt2w  = 0:dt:(ns-1)*dt;               % two way traveltime
%%%%%   Compute t(z) from V(z)
    temptz(1)   = 2.0*firstz/layer_velocity(1);
    for iz=1:nlay-1; 
        temptz(iz+1) = temptz(iz) + ...
            2.0*layer_thickness(iz)/layer_velocity(iz);
    end
    temptz(nlay + 1) = temptz(nlay) + ...
        2.0*(zmax - sum(layer_thickness(1:nlay)))/layer_velocity(nlay);
    tofz  = interp1([0 cumsum(layer_thickness(1:nlay-1)) zmax] ,temptz,z);
%%%%%   Compute z(t) from t(z)
    vfz   = layer_velocity(1);             % velocity of layer 1
    vlz   = layer_velocity(nlay-1);        % velocity of last layer
    lt    = firstt+(ns-1)*dt;
    lz    = firstz+(nz-1)*delz;
    zoft  = yxtoxy(nz,delz,firstz,tofz,ns,dt,firstt,0.0,0.0);
    ii = find(tt2w < tofz(1));             % provide for out of range values
    if ~isempty(ii),
        zoft(ii) = 0.5*tt2w(ii)*vfz;
    end
    ii=find(tt2w >= tofz(nz));
    if ~isempty(ii),
        zoft(ii) = lz + 0.5*(tt2w(ii) - tofz(nz))*vlz;
    end
%%%%%   Compute RMS velocity from z(t) 
    vintt    = 2.0*(zoft(2)-zoft(1))/dt;   % interval velocity
    cum      = firstt*vintt*vintt;
    vrmst(1) = vintt;
    it = 1;   
    t = firstt + dt;
    while it < ns,
        vintt = 2.0*(zoft(it+1)-zoft(it))/dt;
        cum = cum +( dt*vintt*vintt);
        vrmst(it+1) = sqrt(cum/t);
        it = it+1;                                                                    
        t = t+dt;
    end
%%%%%   The Stolt migration velocity is the minimum RMS velocity 
    vmig  = min(vrmst);
%%%%%  Compute stretched-time function
    disp('STOLTMIG > Stretching ... ');
    h = waitbar(0,'Stretching ...');
    scale = 2.0/(vmig*vmig);
    cum   = 0.0;
    st(1)  = cum;
    for it = 1:ns-1,
        cum = cum +( 0.5*dt*(tt2w(it+1)*vrmst(it+1)*vrmst(it+1) + ... 
            tt2w(it)*vrmst(it)*vrmst(it)));
        st(it+1) = sqrt(scale*cum);
    end
%%%%%   Define sampling rate and size of stretched time
    stmax = st(ns);
    delt  = min(diff(st)) / (2.0*fn*dt);
    nss   = 1 + floor(stmax/delt);
%%%%%   Compute t(u) from u(t)
    ts    = yxtoxy(ns-1,dt,0.0,st,nss,delt,0.0,0.0,(ns-1)*dt);
%%%%%   and rid of rerpeating end values in t(u) that will confuse interp1 
    yts=[]; 
    for i=2:nss; 
        if ts(i-1) == ts(i), 
            j=length(ts(1:i-1)); 
            yts=ts(1:j); 
        end; 
    end
    if ~isempty(yts), 
        ts = yts;
        nss = length(yts);
    end
%%%%%    Now time-stretch the section 
    ds = zeros(nss,ntr);
    for i=1:ntr;
        ds(:,i) = interp1(tt2w',d(:,i),ts',interpm,0);
        waitbar(i/ntr, h);
    end
    close(h);
%%%%%   FFT -> FKx
    disp('STOLTMIG > Transforming TX -> FKx');
    fk    = fft2(ds);  
    fk    = fftshift(fk);
%%%%%   Depth sampling interval        
    dz    = zoft(ns)/(nss-1);
end                                        %%%   if nlay loop

%%% Compute image in the time-wavenumber (TKx) domain using the          
img = do_stolt(delt, dx, dz, vmig, fk, nss, ntr, interpm);

%%%%%   Inverse FT KzKx -> ZX and transform -> TX 
disp('STOLTMIG > Transforming KzKx -> ZX -> TX');
ds = ifft2(ifftshift(img));
if nlay == 1,
        dmig = real(ds(1:ns,:));

elseif nlay > 1,
    ds = real(ds(1:nss,:));
    dmig = zeros(ns,ntr);
%%%%%  Unstretch  
    disp('STOLTMIG > Unstretching ... ');
    h = waitbar(0, 'Unstretching ...');
    for i=1:ntr,
        dmig(:,i) = interp1(ts',ds(:,i),tt2w',interpm,0);
        waitbar(i/ntr,h);
    end
    close(h);
end
% Done ...
time  = toc;
disp(['STOLTMIG: finished in ' num2str(time) ' seconds'])
return

function img = do_stolt(dt, dx, dz, vmig, fk, nw, nkx, interpm)
% Computes the image for Stolt F-K migration 

% Angular frequency
w0  = -pi/dt;
dw  = 2*pi/(nw*dt);
w  = [w0 : dw : abs(w0)-dw]';
% Wavenumbers
kx0  = -pi/dx;
dkx  = 2*pi/(nkx*dx);
kx = [kx0 : dkx : -kx0-dkx]';
nkz  = nw;
kz0  = -pi/dz;
dkz  = 2*pi/(nkz*dz);
kz = [kz0 : dkz : -kz0-dkz]';
kzi= [1:1:length(kz)]';

h = waitbar(0,'Computing image FK_x -> K_zK_x');
img = zeros(nw,nkx);
for ikx=1:nkx; 
    ks  = sqrt(kx(ikx)^2 + kz.^2);
    wz  = sign(kz).*ks*vmig /2 ;
    img(:,ikx)=interp1(w, fk(:,ikx), wz, interpm, 0); 
    ij = find(ks);                 % scale the image  
    img(ij,ikx) = img(ij,ikx).*(vmig*abs(kz(ij))./ks(ij));
    ji = setxor(kzi, ij);          % find dc and remove
    img(ji,ikx) = 0;
    waitbar(ikx/nkx, h);
end;
close(h);
return
