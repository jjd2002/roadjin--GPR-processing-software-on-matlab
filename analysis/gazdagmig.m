function dmig  = gazdagmig(d, dt, dx, vofh)
%
% GAZDAGMIG : Gazdag phase-shifting migration for constant or layered
%             velocity structures. 
%       ==>   For the sake of speed, the construction of the image is
%             usually performed with the FORTRAN '90 MEX-file
%             "gazdag.f90". For MS Windows OS, the MEX-file is provided
%             ready to use (gazdag.dll). For Linux OS or any other flavour
%             of Unix, the MEX-file should be built by the user. At any
%             rate,(very much) slower native M-code to perform the same
%             tasks is attached in the subfunctions "gazdagcv" and
%             "gazdaglv". This will take over if the MEX-file is not
%             available.  
%
%    Usage  : dmig  = gazdagmig(d, dt, dx, vofh)
%
%   Inputs  : 
%        d  : The zero- or common-offset GPR section 
%       dt  : The sampling rate 
%       dx  : The trace spacing (spatial sampling rate) 
%  vofh(n,2): The 1-D velocity model of "nlay" velocity -
%             thickness pairs, for example:
%             [ 0.1    1 ;       ... 1st layer
%               0.08   2 ;       ... 2nd layer
%               ...
%               0.18   0 ]       ... n'th layer==basal halfspace.
%              A uniform halfspace is given as a single layer stucture with
%              zero thickness. Velovity values are given in m/ns and
%              thichnesses in m. 
%
%   Outputs  : 
%      dmig  : The migrated section
%
%   Requires : yxtoxy.m
%
%      Uses  : FORTRAN'90 MEX function "gazdag.<mex> with source code
%              "gazdag.f90", and the attached subfunctions "gazdagcv" and
%              "gazdaglv" 
%
%  Author    : Andreas Tzanis,
%              Dept. of Geophysics,
%              University of Athens
%              atzanis@geol.uoa.gr
%
%  Copyright (C) 2005, Andreas Tzanis. All rights reserved.
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


[ns,ntr]=size(d);
%%%%%   Check Velocity structure
[nlay, vhpairs] = size(vofh);
if vhpairs ~=2,
    erh = errordlg('Error in the structure of the velocity model',...
        'GAZDAGMIG : ERROR');
    uiwait(erh)
    dmig = [];
    return
end
layer_velocity  = vofh(:,1)';
layer_thickness = vofh(:,2)';
for i=1:nlay
    if layer_velocity(i) <= 0 || layer_velocity(i) > 0.2998,
        errordlg('Impossible velocity value found! Please try again!',...
            'GAZDAGMIG : ERROR')
        dmig = [];   
        return
    end
end
%%%  Begin
tic                                      % start counting processing time
if nlay == 1,                            %%%%%   Case of uniform halfspace
    vmig = layer_velocity(1);
elseif nlay > 1,                         %%%%%   Case of layered halfspace
    tt2w      = 0:dt:(ns-1)*dt;          % two way traveltime
    firstz = 0.0;
    firstt = 0.0;
%%% Compute migration velocity vmig(t) from v(z) 
    zmax = max(vofh(:,1))*tt2w(ns)/2;    % max possible penetration
    delz = (zmax-firstz)/(ns-1);
    z    = firstz : delz : zmax;
    nz        = length(z);
%%% Compute t(z) from V(z)
    temptz(1) = 2.0*firstz/layer_velocity(1);
    for iz=1:nlay-1; 
        temptz(iz+1) = temptz(iz) + ...
            2.0*layer_thickness(iz)/layer_velocity(iz);
    end
    temptz(nlay + 1) = temptz(nlay) + ...
        2.0*(zmax - sum(layer_thickness(1:nlay)))/layer_velocity(nlay);
    tofz  = interp1([0 cumsum(layer_thickness(1:nlay-1)) zmax] ,temptz,z);
%%% Compute z(t) from t(z)
    vfz   = layer_velocity(1);             % initial velocity
    vlz   = layer_velocity(nlay-1);        % final velocity at depth z
    lt    = firstt+(ns-1)*dt;
    lz    = firstz+(nz-1)*delz;
    zoft  = yxtoxy(nz,delz,firstz,tofz,ns,dt,firstt,0.0,0.0);
    ii = find(tt2w < tofz(1));             % out of range values
    if ~isempty(ii),
        zoft(ii) = 0.5*tt2w(ii)*vfz;
    end
    ii=find(tt2w >= tofz(nz));
    if ~isempty(ii),
        zoft(ii) = lz + 0.5*(tt2w(ii) - tofz(nz))*vlz;
    end
%%% Compute vmig(t) from z(t) 
    vmig = [];
    for it = 1:ns-1,
        vmig(it) = 2.0*(zoft(it+1) - zoft(it))/dt;
    end
    vmig(ns) = vmig(ns-1);
    disp('GAZDAGMIG > Computed velocity, OK');
end                                         %%%%%   if nlay loop

%%% Fourier transform to frequency-wavenumber (FKx) domain
fk  = fft2(d);      
fk  = fftshift(fk);
try 
% Compute the image in using the fortran MEX  routine
    disp('GAZDAGMIG > Calling MEX-file "gazdag"');
    if nlay == 1,
        [imgr, imgi] = gazdag(dt, dx, vmig, real(fk), imag(fk), ...
            ns, ntr, 'uniformv');  
    elseif nlay >1 
        [imgr, imgi] = gazdag(dt, dx, vmig, real(fk), imag(fk), ... 
            ns, ntr, 'layeredv');
    end                            
    img = imgr +sqrt(-1)*imgi;      %%% Get complex image
    clear imgr imgi;               
%%% Alternatively, use the slower MATLAB code in the attached subfunctions 
%%% gazdagcv and gazdaglv, but only if the MEX file is invalid (unsupported
%%% compilers) or missing 
catch ME %#ok<NASGU>
    idSegLast = regexp(ME.identifier, '(?<=:)\w+$', 'match');
    if strcmp(idSegLast,'invalidMEXFile')
        disp('GAZDAGMIG > MEX-file "gazdag" invalid or not found!');
        disp('GAZDAGMIG > Reverting to slower MATLAB code!');
        if nlay == 1,
            img = do_gazdagcv(ns, ntr, dt, dx, vmig, fk);
        elseif nlay >1
            img = do_gazdaglv(ns, ntr, dt, dx, vmig, fk);
        end
%%% Otherwise, terminate program and notify user
    else 
        rethrow(ME); 
    end 
end
%%% Transform from time-wavenumber (TKx) to time-space (TX) domain and 
%%% get migrated section
h = waitbar(0,'Transforming TK_x -> TX');
dmig = zeros(ns,ntr);
for i=1:ns,
    dmig(i,:) = ifft(ifftshift(img(i,:)));
    waitbar(i/ns,h);
end
dmig = real(dmig);
close(h);
time  = toc;
disp(['GAZDAGMIG finished in ' num2str(time) ' seconds'])
return
%
function img = do_gazdagcv(ns, ntr, dt, dx, vmigc, fk)
%
% This is the MATLAB code to compute the FKx -> TKx image of Gazdag
% Phase-shifting migration for zero-offset GPR data in uniform halfspaces
% (migration velocity  vmigc is constant).
% *** The active implementation is vectorized (as much as the
%     integration process will allow and is reasonably fast. It will do a
%     [512 x 512] job in about 36 seconds (as opposed to the 5 seconds
%     needed for the MEX-file which works with typical -FORTRAN-
%     sequential programing). 
% *** The inactive (commented out) implementation is traditional sequential
%     programing (the same as in FORTRAN). It is VERY SLOW!
% *** Note that the continuation operator "cc" is conjugated, opposite to
%     what books usually say), to account for the (engineering) definition
%     of the Fourier kernel in MATLAB  
%
% Author: Andreas Tzanis (C) 
%         November 2005
%
% Set up physical constants 
nw   = ns;     w0   = -pi/dt;     dw   = 2*pi/(ns*dt);
nx   = ntr;    kx0  = -pi/dx;     dkx  = 2*pi/(nx*dx);
ntau = ns;     dtau = dt;
% initialize the image array
img = zeros(ns,ntr); 
% Preare the progress bar
h=waitbar(0,'Imaging FK_x -> TK_x');
% The following two lines needed for the vectorized program ---------------
kx = kx0 : dkx : -kx0-dkx;
kx2 = kx.*kx;
% -------------------------------------------------------------------------
for iw = 1:nw,
    w  = w0 + (iw-1)*dw;
    if w == 0.0,
        w = 1e-10/dt;
    end
	progress = iw/nw;                    % Report progress 
    waitbar(progress,h);
% THIS IS THE SEQUENTIAL PROGRAMMING AS IN FORTRAN ------------------------
%    for ikx = 1:nx,
%        kx  = kx0 + (ikx-1)*dkx;
%        w2   = w*w;
%        vkx2 = (vmigc*vmigc * kx*kx)/4.;
%        if w2 > vkx2,
%            phase = real(-w * dtau * sqrt(1.0d0 - vkx2/w2));
%            cc    = conj(complex(cos(phase),sin(phase)));
%            % Accumulate image summed over all frequencies
%            for itau = 1:ntau,
%                fk(iw,ikx)   = fk(iw,ikx) * cc;
%                img(itau,ikx)= img(itau,ikx) + fk(iw,ikx);
%            end
%        else
%            fk(iw,ikx) = complex(0.0,0.0);
%        end
%    end                                            % ikx loop
% THIS IS THE SEMI-VECTORIZED PROGRAMMING ---------------------------------
    w2   = w*w;
    vkx2 = (vmigc*vmigc * kx2)/4.;
    ik = find (vkx2 < w2);
    ffk = fk(iw,ik);
    phase = real(-w * dtau * sqrt(1.0d0 - vkx2(ik)/w2));
    cc    = conj(complex(cos(phase),sin(phase)));
    % Accumulate image summed over all frequencies
    for itau = 1:ntau,
         ffk = ffk.*cc;
         img(itau,ik)= img(itau,ik) + ffk;
    end
end                                                % iw loop
% Normalize for inverse FFT
img = img/nw;
% Close the progress bar
close(h); 
return
%
function img = do_gazdaglv(ns, ntr, dt, dx, vmigv, fk)
%
% This is the MATLAB code to compute the FKx -> TKx image of Gazdag
% Phase-shifting migration for zero-offset GPR data in LAYERED halfspaces
% (migration velocity  vmigc is a function of time).
% *** The active implementation is (to a point) vectorized and executes
%     reasonably fast. It will do a [512 x 512] job in about 97 seconds. 
% *** This routine is programmed sequentially, as in FORTRAN (same code
%     with gazdag.f90 and is DAMNED SLOW!
% *** Note that the continuation operator "cc" is conjugated, (opposite to
%     what books usually say), to account for the (engineering) definition
%     of the Fourier kernel in MATLAB   
%
% Author: Andreas Tzanis (C) 
%         November 2005
%

nw   = ns;      w0   = -pi/dt;      dw   = 2*pi/(ns*dt);
nx   = ntr;     kx0  = -pi/dx;      dkx  = 2*pi/(nx*dx);
ntau = ns;      dtau = dt;
ft   = 0;   	ftau = ft;      	tmax  = ft  + (ntau-1)*dtau;
img = zeros(ns,ntr); 
% Preare the progress bar
h=waitbar(0,'Imaging FK_x -> TK_x');
% THIS IS THE SEQUENTIAL PROGRAMMING AS IN FORTRAN ------------------------
% for ikx = 1:nx                                     % Loop over wavenumbers
%     kx  = kx0 + (ikx-1)*dkx;
% 	prgrs = ikx/nx;                                 % Report progress  
%     waitbar(prgrs,h);
%     for itau = 1:ntau,                             % Loop over time
%          tau = ftau + (itau-1)*dtau;
%          for iw = 1:nw,                            % Loop over frequencies
%              w  = w0 + (iw-1)*dw;
%              if w == 0.0,
%                  w = 1e-10/dt;
%              end
%              coss = 1.0 - (0.5 * vmigv(itau)*kx/w)^2;
%              if coss > (tau/tmax)^2,
%                  phase = real(-w*dt*sqrt(coss));
%                  cc    = conj(complex(cos(phase),sin(phase)));
%                  fk(iw,ikx)   = fk(iw,ikx) * cc;
%              else
%                  fk(iw,ikx) = complex(0.0,0.0);
%              end   
%              img(itau,ikx)= img(itau,ikx) + fk(iw,ikx);
%          end                                     % iw loop
%      end                                         % itau loop
%  end                                    
% THIS IS THE SEMI-VECTORIZED PROGRAMMING ---------------------------------
w = w0 : dw : -w0-dw;
iw = find(w==0);
w(iw) = 1e-10/dt; %#ok<FNDSB>
w = w(:);
for ikx = 1:nx                                      % Loop over wavenumbers
    kx  = kx0 + (ikx-1)*dkx;
	prgrs = ikx/nx;                                       % Report progress  
    waitbar(prgrs,h);
    for itau = 1:ntau,                                     % Loop over time
         tau = ftau + (itau-1)*dtau;
         coss = 1.0 - (0.5 * vmigv(itau)*kx./w).^2;
         iw = find((tau/tmax)^2 < coss);
         phase = real(-w(iw).*dt.*sqrt(coss(iw)));
         cc    = conj(complex(cos(phase),sin(phase)));
         fk(iw,ikx)   = fk(iw,ikx) .* cc;
         iw = find((tau/tmax)^2 >= coss);
         fk(iw,ikx) = complex(0.0,0.0); %#ok<FNDSB>
         img(itau,ikx)= img(itau,ikx) + sum(fk(:,ikx));
    end                                                         % itau loop
end
img = img/nw;
% Close the progress bar
close(h);
return
