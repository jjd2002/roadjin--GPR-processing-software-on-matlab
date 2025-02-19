function tx = inv_taup_tx(tp, dt, x, dx, p)
% 
% INV_TAUP_TX : Compute an inverse Tau-P transform (inverse slant stack)
% 				in T-X domain
%
% Usage : tx = inv_taup_tx(tp, dt, x, dx, p);
%
%  Input:
%    tp : 2-D array of input traces in Tau-P domain
%    dt : time sampling interval
%     x : vector of receiver (trace) locations
%    dx : trace spacing in x
%     p : vector of slopes (slowness) for Tau-P transform
% 
% Output:
%    tx : 2-D array of output traces in T-X domain
%
% Credit: Programming tips from C old code written by Gabriel Alvarez, CWP
%
% Author: Andreas Tzanis,
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
% *************************************************************************

% Useful parameters
nt    = size(tp,1);
np    = max(size(p));
nx    = max(size(x));

% Design and apply rho-filter
%%%%%%%%%%%%    Time domain filtering    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho = rho_filter(127,nt,dt);
tpr = convn(tp,rho','same');
%%%%%%%%%%%%    Frequency domain filtering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nfft = 2*(2^nextpow2(nt));
% nfh  = nfft/2;  
% df   = 1.0/(2*dt*(nfh));
% cf   = [ complex(df*(1:nfh), 0.0) fliplr(complex(df*(1:nfh), 0.0))];
% TPR  = fft(tp,nfft,1);
% TPR  = TPR .* (cf' * ones(1, np));
% tpr  = real(ifft(TPR,nfft,1));
% tpr  = tpr(1:nt,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Invert Tau-P to T-X
tx = zeros(nt, nx);                                % Initialize ouput array
p  = fliplr(p);
barhandle = waitbar(0,'Transforming Tau-P --> T-X');
for ix = 1:nx,
    for ip = 1:np
        delay = p(ip)*x(ix)/dt;
        if delay >=0,
            id = fix(delay);
            i1 = id+1;
            i2 = nt - 1;
        else
            id = fix(delay)-1;
            i1 = 1;
            i2 = nt+id;
        end
        fraction = delay-id;
        for it = i1:i2-1
            tx(it-id,ix) = tx(it-id,ix) + ...
                abs(dx)*(tpr(it,ip) + fraction*(tpr(it+1,ip) - tpr(it,ip)));
        end
    end
    waitbar(ix/nx,barhandle);
end
close(barhandle)
return

function [rho, cf] = rho_filter(nf,nt,dt)
% 
% RHO_FILTER : Compute the rho filter in frequenccy domain for the time
%              domain inverse Tau-P transform
% 
% Input:
%   nf : number of point for the rho filter
%   nt : number of time samples
%   dt : time sampling interval
% 
%Output:
%  rho : 1-D array of filter points
% 
% Author: Andreas Tzanis,
%         Department of Geophysics, 
%         University of Athens
%         atzanis@geol.uoa.gr
%
% Copyright (C) 2008, Andreas Tzanis. All rights reserved.
%
% *************************************************************************

% Ensure that number of filter points is odd
if rem(nf,2) == 0,
    error('==> RHO_FILTER > Filter length mu be odd number!');
end
nph = round(nf/2);
% Compare filter length with number of time samples 
if nph > nt,
    error('==> RHO_FILTER > filter length larger than number of time samples!');
end
nfft = 2*(2^nextpow2(nt));                       % compute padding factor 
% compute filter coefficients 
nfh = nfft/2;
df = 1.0/(2*dt*nfh);
cf = [ complex(df*(1:nfh), 0.0) fliplr(complex(df*(1:nfh), 0.0))];
% inverse Fourier transform from f to t 
cf  = ifft(cf);
rho = zeros(1,nf);
% normalize output filter coefficients 
rho(nph) = 1.0;
for it = 1:nph-1,
    rho(nph - it) = real(cf(it + 1))/real(cf(1));
    rho(nph + it) = real(cf(it + 1))/real(cf(1));
end
return
