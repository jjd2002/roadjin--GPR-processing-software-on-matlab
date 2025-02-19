function [yo, blx, bly] = baseline(istep,y,action,interpm);
%
% BASELINE : Determine a baseline BLY(NBL) for the 
%            vector Y(N) and optionally remove it. 
%
%    Usage : [yo, blx, bly] = baseline(istep,y,action,interpm)
%
%   Inputs : 
%    istep  = bin size for determining the baseline 
%        y  = the data vector
%  action   = 'base' determine a baseline ONLY
%           = 'high' determine and remove a baseline (high-pass filter)
%           = 'low'  determine baseline and replace the input series
%             with the baseline(low pass filter)
% interpm     defines the interpolation method to be used. Any of the
%             methods implemented in MATLAB's "interp1" routine can be
%             used, i.e. 'nearest', 'linear', 'spline', 'pchip', 'cubic' 
%             and 'v5cubic'. Default is linear interpolation.
%
%  Outputs :    
%       yo  = the processed data vector 
%      blx  = x-coordinates of the base line w.r.t. the elements of the
%             input vector
%      bly  = the baseline in input data units
%
%  Author : Andreas Tzanis,
%           Department of Geophysics, 
%           University of Athens
%           atzanis@geol.uoa.gr
%
% (C) 2005, Andreas Tzanis, all rights reserved
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
if nargin < 4, 
    interpm = 'linear';       % default interpolation method 
end
interpm = lower(interpm);
if nargin < 3, 
    action  = 'base';         % default is to determine baseline only
end
action  = lower(action);

n   = length(y);
yo  = zeros(n,1);
nsteps = floor(n/istep); 
%% Determine the baseline
blx = [1];
bly = [y(1)];
for i=1:nsteps
    ii = (i-1)*istep +1;
    blx = [blx; ii+floor(istep/2) ];
%    bly = [bly; mean(y(ii : ii+istep-1))];
    bly = [bly; median(y(ii : ii+istep-1))];
end
ii = n - istep*nsteps;
if ii == 0, 
    blx = [blx; n];
    bly = [bly; y(n)];
else
    blx = [blx; n];
%    bly = [bly; mean(y(n-ii+1 : n))];
    bly = [bly; median(y(n-ii+1 : n))];
end

% High pass filter mode - Remove the baseline
if strcmp(action,'high'),
    nbl = length(blx);
    yo  = interp1(blx, bly, [blx(1) : 1 : blx(nbl)]', interpm, 0);
    yo = y - yo;
    return;
end

%%% Low-pass filter mode - interpolate the baseline
if strcmp(action,'low')
    nbl = length(blx);
    yo  = interp1(blx, bly, [blx(1) : 1 : blx(nbl)]', interpm, 0);
    return
end
return;
