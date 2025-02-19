function dm = mat_radar2d(K,M,Q,nx,nz,dx,dz,nt,dt,freqc,band)
%
%*******************************************************************
%  SPLIT-STEP 2D GPR MODELLING after 
%  Grandjean, G., 1998, Geophysical Prospecting, 46, 287-301, and,
%  Grandjean, G. and Durand, H., 1999,Computers and Geosciences, 25 141-149.                   *
%********************************************************************
%
% Imports matrices of model relative dielectric constant, magetic
% permeability and quality factor structures and computes the corresponding
% velocity structure after Bano (1996, Geophys. J. Int., 279 - 288) and
% hence the forward model using a conjugate split-step method. 
% ==> This program is part of the MATGPR suite.
%
% Usage  : dm = mat_radar2d(K,M,Q,nx,nz,nt,dx,dz,dt,freqc,band);
%
%  Input : K            : Matrix of the model relative dielectric constant
%                         with dimensions (nz x nx)
%          M            : Matrix of the model magnetic permeability with
%                         dimensions (nz x nx)
%          Q            : Matrix of the model quality factor with
%                         dimensions (nz x nx)
%          nx           : Model dimension along the scan axis
%          nz           : Vertical model dimension (along depth axis)
%          dx           : Trace spacing of the synthetic model
%          dz           : Depth spacing of the synthetic model
%          nt           : Dimension of the time axis in the synthetic
%                         radargram
%          dt           : The sampling rate in the synthetic radargram 
%          freqc        : The central frequency of the GPR antenna.
%          band         : Bandwidth selection flag
%                         = 0, limited frequency scan, 0.25*freqc<f<2*freqc
%                         ~=0, full frequency scan, 0 < f <= Nyquist 
%
% Output : dm           : The synthetic GPR section 
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

%%%% Determine frequency range  
m0   = 4*pi*1.0e-7;
e0   = 8.8592e-12;
nx2    = nx;
nt2    = 2^nextpow2(nt);
mhz    = 1000000;
tet1   = 60;
tet2   = 80;
wr     = 2.0*pi*freqc*mhz;
freq1  = fix(freqc*mhz*0.25);
freq2  = fix(freqc*mhz*2);
freqc2 = fix(freqc*mhz);
iw1    = fix(nt2*dt*freq1);
iw2    = fix(nt2*dt*freq2);
 if band < 0.0 || band > 0.0, 
    iw1 = 1;
    iw2 = (nt2/2) +1;
 end
iwc    = fix(nt2*dt*freqc2);
%report status
disp('> Computed parameters, OK.');

% Initialize arrays
zp     = zeros(nz,nx2);
cp     = zeros(nt2,nx2);
dm     = zeros(nt,nx);
waven  = zeros(1,nx2);
omega  = zeros(iw2,1);
g      = zeros(iw2,1);

%%%% COMPUTE VELOCITY STRUCTURE 
vit = 1.0./(sqrt(m0.*M*e0.*K).*cos((pi/4).*(1.0-((2/pi).*atan(Q)))));
%%%% REFLECTION COEFFICIENTS
for ix = 1:nx,
   for iz = 2:nz
         cref1 = sqrt( complex(M(iz,ix), 0.0) / ...
                 complex( K(iz,ix)*sin(atan(Q(iz,ix))), ...
                 cos(K(iz,ix)*atan(Q(iz,ix)))));
         cref2 = sqrt( complex(M(iz-1,ix),0.0) / ...
                 complex( K(iz-1,ix)*sin(atan(Q(iz-1,ix))), ...
                 cos(K(iz-1,ix)*atan(Q(iz-1,ix)))));
         zp(iz,ix) = (cref1 - cref2)/(cref1 + cref2);
   end
end
%%%% mean velocity 
vitm = mean(vit/2,2);
%%%% mean quality factor
qm   = mean(Q,2);
disp('> Computed velocity structure, OK.');                 %report status

%%%% Angular frequencies and wavenumbers
for ikx = 1:nx2
   kx = ikx*2.0*pi/nx2;
   if kx > pi,
       kx = 2.0*pi - kx;
   end
   waven(ikx) = kx/dx;
end
for iw = iw1:iw2,
   omega(iw) = (iw*2.0*pi/nt2)/dt;
end 
if band > 0.0,
    sigma = (iw2-iw1)*band;
    for iw = iw1:iw2,
        g(iw) = exp(-1.0*((iw-iwc)/sigma)*((iw-iwc)/sigma));
    end
end
% Topmost wavenumber^2 for radiation pattern calculations
bomega2 = (omega/vitm(1)).*(omega/vitm(1));
disp('> Computed frequencies and wavenumbers, OK.');         %report status

%%%% BEGIN FORWARD COMPUTTIONS
%%%% Transform X -> Kx
zp = fft(zp,[],2);
disp('> Forward computations: Initial FT X -> Kx, OK.');
hb = waitbar(0,'Computing model ... please wait!');
%%%%  LOOP OVER DEPTHS
for iz = nz:-1:1,
    waitbar((nz-iz)/nz,hb);                               % Report progress
    nq = (2.0/pi)*atan(qm(iz));
    expq = (1.0 - nq)/2.0;
    % Assign propagation constants for this depth 
    bomega = omega./ (vitm(iz) * (omega/wr).^expq);
    aomega = bomega*tan((1.0 - nq)*(pi/4.0));
    komega = complex(bomega, aomega);
    denc   = komega.*komega;
    %%%%  LOOP OVER WAVENUMBERS AND FREQUENCIES
    for ikx = 1:nx2,
        den = waven(ikx)*waven(ikx);
        denkc = complex(den, 0.0);
        denkz = sqrt(denc - denkc);
        %%%%  INCLUDE RADIATION PATTERN
        difc  = bomega2 - den; 
        kz    = sqrt(difc);
        tet   = atan(waven(ikx)./kz);
        radc  = ones(iw2,1);
        ii    = find(difc >0);
        if ~isempty(ii),
            alf1 = pi*tet1/180.0;
            alf2 = pi*tet2/180.0;
            jj   = find(alf1 < tet(ii) & tet(ii) <= alf2);
            if ~isempty(jj),
                radc(ii(jj)) = ...
                    (0.42 - 0.5*cos((1.0 + ...
                    (tet(ii(jj))-alf1)/(alf2-alf1))*pi)+  ...
                    0.08*cos((1.0 + ...
                    (tet(ii(jj))-alf1)/(alf2-alf1))*2.0*pi));
            end
            jj   = find(alf2 < tet(ii) );
            if ~isempty(jj),
                radc(ii(jj)) = 0.0;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        den1 = real(denkz);
        den2 = imag(denkz);
        rshift = exp(1.0*dz*den2);
        ii = find(den2 > 0 );
        if ~isempty(ii)
            rshift(ii) = exp(-1.0*dz*den2(ii));
        end
        shift1 = exp(complex(0.0, dz*den1));
        shift2 = shift1.*rshift;
        cp(iw1:iw2,ikx) = cp(iw1:iw2,ikx).*shift2;
        if band > 0.0,
            cp(iw1:iw2,ikx) = cp(iw1:iw2,ikx).*g;
        end
        cp(iw1:iw2,ikx) = cp(iw1:iw2,ikx) + zp(iz,ikx);
        cp(iw1:iw2,ikx) = cp(iw1:iw2,ikx).*radc;
    end
    %%%% Inverse transform Kx -> X
    cp(iw1:iw2,1:nx2) = ifft(cp(iw1:iw2,1:nx2),[],2);
    %%%%   CORRECTION FOR LATERAL VELOVITY VARIATIONS
    for ix = 1:nx,
        nq = (2.0/pi)*atan(Q(iz,ix));
        expq = (1.0 - nq)/2.0;
        bomegaxz = omega./ ((vit(iz,ix)*((omega/wr).^expq) )/2.0);
        aomegaxz = bomegaxz*tan((1.0-nq)*(pi/4.0));
        komegaxz = complex(bomegaxz,aomegaxz);
        phase = complex(0.0,-1.0)*(komega - komegaxz);
        shift = exp(phase*dz);
        cp(iw1:iw2,ix) = cp(iw1:iw2,ix).*shift;
    end
    %%%%  Transform X -> Kx
    cp(iw1:iw2,1:nx2) = fft(cp(iw1:iw2,1:nx2),[],2);
end                                                %%% END LOOP OVER DEPTHS
close(hb)
disp('> Finished forward computations, OK.');
%%%% Inverse transform Kx -> X
cp(iw1:iw2,1:nx2) = ifft(cp(iw1:iw2,1:nx2),[],2);
disp('> Finished inverse FT Kx -> X');
%%% DONE - transform from frequency to time domain
for ix = 1:nx,
    for iw = iw1:iw2,
        cp(nt2-iw+1,ix) = conj(cp(iw,ix));
    end
    cp(:,ix) = fft(cp(:,ix),nt2);
    dm(1:nt,ix) = real(cp(1:nt,ix));
end
disp('> Done inverse FT K -> T OK - OVER and OUT.');
return
