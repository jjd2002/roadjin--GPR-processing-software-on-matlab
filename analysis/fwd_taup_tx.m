function tp = fwd_taup_tx(d, dt, x, dx, p)
% 
% FWD_TAUP_TX : Compute a forward taup transform (slant stack) in the 
% 			    T-X domain
%
%  Usage: tp = fwd_taup_tx(d, dt, x, dx, p);
%
%  Input:
%     d : 2-D array of input traces in T-x domain
%    dt : time sampling interval
%     x : vector of receiver (trace) locations
%    dx : trace spacing in x
%     p : vector of slopes (slowness) for Tau-P transform
% 
% Output:
%    tp : 2-D array of output traces in Tau-P domain
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

[nt, nx] = size(d);
np       = max(size(p));
tp       = zeros(nt,np);                           % Initialize ouput array
barhandle = waitbar(0,'Transforming T-X -> Tau-P');
for ip = 1:np,
    for ix = 1:nx
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
            tp(it-id,ip) = tp(it-id,ip) + ...
                abs(dx)*(d(it,ix) + fraction*(d(it+1,ix) - d(it,ix)));
        end
    end
    waitbar(ip/np,barhandle);
end
close(barhandle)
return
