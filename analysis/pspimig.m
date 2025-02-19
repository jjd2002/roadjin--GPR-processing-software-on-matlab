function [dmig, zmig, nz] = pspimig(d, dt, dx, v2d, zv, dz)
%
% PSPIMIG : Migration of zero-offset data by phase shift plus
%           interpolation after Gazdag, J. and Sguazzero, P., 1984.,
%           Migration of seismic data by phase shift plus interpolation,
%           Geophysics 49(2), 124-131.
%
%   Usage : [dmig, zmig, nz] = pspimig(d, dt, dx, v2d, zv, dz)
%
%  Inputs :
%       d : The [ns x ntr] array of the equally-spaced GPR section
%      dt : The sampling rate in the ns dimension (in nanosecs)
%      dx : Trace spacing in the ntr dimension (in m).
%     v2d : The velocity model. 
%           ** If v2d is a single matrix, then it must be the [nv x ntr]
%              non-dispersive term of the velocity model. 
%           ** If v2d is a {1 x 1} or a {1 x 2} cell array, then v2d{1}
%              must be the [nv x ntr] non-dispersive term of the velocity
%              model, which is the only term used by PSPIMIG. 
%      zv : [nv x 1] vector, the depth coordinates of the velocity model.
%      dz : Depth spacing used for integration (stepping size).
%
% Outputs : 
%    dmig : [nz x ntr] array, the depth migrated GPR section.
%    zmig : [nz x 1] vector, the depth coordinates of the migrated section.
%      nz : The vertical dimension of dmig == dimension of zmig
%
%Requires : checkcomma.m
%
%  Author : Andreas Tzanis, 
%           Department of Geophysics, 
%           University of Athens, 
%           atzanis@geol.uoa.gr
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

if nargin < 6,           % Depth stepping size may not have been furnished
    dz = 0.02; 
    cldef=cellstr(num2str(dz));
    answer = inputdlg('Give depth spacing in m', ...
        'MIGPSPI : Requesting dz',1 , cldef); 
    if isempty(answer)
        dmig = [];   zmig = [];    nz = [];
        return      
    end
    lb = checkcomma(answer);     
    dz = str2num(lb(1,:));
end
%%% Velocity model and its depth coordinates not supplied
if nargin < 5 || nargin < 4,
    info = {'Velocity information not supplied or incomplete.' ...
          'Please generate a 2D velocity model using GET2DVELOCITY'}; 
    erh = errordlg( info,'MIGPSPI: ERROR');
    uiwait(erh)
    dmig = [];    zmig = [];    nz=[];
    return
end

% Extract the velocity model from input cell array.
if ~iscell(v2d),
    v = v2d;
else
    v = v2d{1};
end
clear v2d;

[nt,nx] = size(d);
% Trap an obvious error ...
if nx ~= size(v,2), 
    erh = errordlg('Mismatch between the sizes of data and velocity model',...
        'MIGPSPI: ERROR');
    uiwait(erh)
    dmig = [];    zmig = [];    nz = [];
    return
end

% Now proceed - Get zmax from user
zmax = zv(max(size(zv))); 
cldef(1) = cellstr(num2str(zmax));
answer = inputdlg('Give maximum depth to step to, in m',...
    'MIGPSPI : INPUT ',1 , cldef); 
if isempty(answer)
    dmig = [];   zmig = [];    nz = [];
    return      
end
lb = checkcomma(answer);    
zmax = str2num(lb(1,:));
%%%% Depth vector
zmig    = 0 : dz : zmax;
nz      = length(zmig);                            % number of depth steps

% determine frequency sampling for t -> w FFT
ntfft = nt;    %ntfft = 2^nextpow2(nt);
nw = floor(ntfft/2) + 1;
dw = 2.0*pi/(ntfft*dt);
fw = 0;
w  = (fw : dw : (nw-1)*dw);
w(1) = 1.0e-10/dt;
% determine wavenumber sampling for x -> k FFT
nxfft = 2^nextpow2(nx);
nk = nxfft;
% dk = (1/(2*dx))/((nk/2)-1); 
% fk = -(1/(2*dx)); 
% kx = 2*pi*(fk : dk : (nk/2)*dk); 
dk = 2.0*pi/(nxfft*dx);
fk = -pi/dx;
kx = fk : dk : (nk/2 -1)*dk;
% Wavenumber and frequency matrices for vector processing
[K, W] = meshgrid(kx, w);
clear w kx;

% allocate space for results
dmig = zeros(nz,nx);
%%%  Begin
tic                                        % start counting processing time
% Determine reference velocities
vmax   = max(max(v));
vmin   = min(min(v));
rhomax = 1.5;
R      = vmax/vmin;
if R >= rhomax,
    L    = ceil(log(R)/log(rhomax) +1);
    rho  = R^(1/(L-1));
    V    = zeros(L,1);
    V(1) = vmin;
    for l=2:L-1
        V(l) = V(l-1)*rho;
    end
else
    L = 2;
    V = [vmin; vmax];
end
V(L) = vmax*1.00005;  

% Transform t --> w
cd = fft(d,ntfft,1);
cd = cd(1:nw,:);

%loops over depth
hwbar   = waitbar(0,['Imaging FK_x -> ZX in ' int2str(nz) ' depth steps']);
for iz = 1:nz,
    % get velocity at iz'th depth
    izv = find(zv <= zmig(iz));
    izv = izv(length(izv));
    
    % Accumulate the migrated data
    for ix = 1:nx,
        for iw = 1:nw,
            dmig(iz,ix) = dmig(iz,ix) + real(cd(iw,ix));
        end;
    end
    
    % Velocity extrema at this depth ...
    lvmax = max(v(izv,:));
    lvmin = min(v(izv,:));

    % Initial phase shift
    for ik = 1:nx,
        cd(1:nw,ik) = cd(1:nw,ik).*exp(-sqrt(-1)*W(:,1)*dz*2.0/v(izv,ik));
    end
    % Transform W-x to W-K domain
    aux = fftshift(fft(cd,nk,2),2);
    
    % The second time phase shift
    v1=lvmin*0.5;
    v2=lvmax*0.5;
    if (v2-v1)/v1 < 0.01,                    % case of small velocity range

        kz = 1.0 - (v1*K./W).^2;
        phase   = zeros(nw,nk);
        cshift  = zeros(nw,nk);
        kplus   = find(kz > 0);
        phase(kplus)  = -W(kplus).*sqrt(kz(kplus))*dz/v1 + W(kplus)*dz/v1;
        cshift(kplus) = complex(cos(phase(kplus)), sin(phase(kplus)));
        aux(kplus) = aux(kplus).*cshift(kplus);
        kminus  = find(kz < 0);
        phase(kminus)  = -W(kminus).*sqrt(-kz(kminus))*dz/v1;
        cshift(kminus) = exp(complex(phase(kminus), W(kminus)*dz/v1));
        aux(kminus) = aux(kminus).*cshift(kminus);
        % Inverse transform W-K to W-x
        cd = ifft(ifftshift(aux,2),nk,2);
        cd = cd(1:nw,1:nx);

    else                               % case of significant velocity range

        aux2 = zeros(nw, nk, L);
        for il = 1:L,                      % Loop over reference velocities
            kz = 1.0 - ((V(il)/2.0)*K./W).^2;
            phase   = zeros(nw,nk);
            cshift  = zeros(nw,nk);
            aux1     = zeros(nw,nk);
            kplus   = find(kz > 0);
            phase(kplus)  = -W(kplus).*sqrt(kz(kplus))*dz*2.0/V(il) ...
                + W(kplus)*dz*2.0/V(il);
            cshift(kplus) = complex(cos(phase(kplus)), sin(phase(kplus)));
            aux1(kplus) = aux(kplus).*cshift(kplus);
            kminus  = find(kz < 0);
            phase(kminus)  = -W(kminus).*sqrt(-kz(kminus))*dz*2.0/V(il);
            cshift(kminus) = exp(complex(phase(kminus), W(kminus)*dz*2.0/V(il)));
            aux1(kminus) = aux(kminus).*cshift(kminus);
            aux2(:,:,il) = aux1;
            % Inverse transform W-K to W-x
            aux2(:,:,il) = ifft(ifftshift(aux2(:,:,il),2),nk,2);
        end

        % Interpolation ...
        for ix = 1:nx,
            for il = 1:L-1,
                if ((v(izv,ix) >= V(il)) && (v(izv,ix) < V(il+1))),
                    v1 = V(il);
                    v2 = V(il+1);
                    a1 = real(aux2(:,ix,il));
                    a2 = real(aux2(:,ix,il+1));
                    b1 = imag(aux2(:,ix,il));
                    b2 = imag(aux2(:,ix,il+1));
                    a = a1*(v2-v(izv,ix))/(v2-v1) + a2*(v(izv,ix)-v1)/(v2-v1);
                    b = b1*(v2-v(izv,ix))/(v2-v1) + b2*(v(izv,ix)-v1)/(v2-v1);
                    cd(:,ix) = complex(a,b);
                    break;
                end
            end
        end
    end                       % END LOOP: if (v2-v1)/v1 < 0.01 ... else ...
    waitbar(iz/nz, hwbar);
end                                                   %END LOOP OVER DEPTHS
dmig = dmig/ntfft;

time  = toc;
disp(['MIGPSPI finished in ' num2str(time) ' seconds'])
close(hwbar)
return;
