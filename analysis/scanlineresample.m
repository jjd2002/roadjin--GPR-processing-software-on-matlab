function [dout, dxout, ntrout, xout] = scanlineresample(d, dx, x)
%
% SCANLINERESAMPLE : Driver for routine RESAMPLE1.M - changes the spatial
% samling rate of the GPR section d to a shorter or longer trace spacing. 
%
%   Usage : [dout, dxout, ntrout] = scanlineresample(d, dx)
%
%  Inputs : 
%     d   : Input 2-D GPR data
%     dx  : Trace spacing in m
%      x  : Trace locations along the scan line
% 
% Outputs :   
%   dout  : The resampled GPR data
%   dxout : The update sampling rate
%    xout : The updated trace location vector 
%  ntrout : The number of traces dout
%
%REQUIRES : resample1.m
%
%  Author : Andreas Tzanis, 
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

if nargin < 3 || ~exist('x') || isempty(x),
    x0 = 0;
else
    x0 = x(1);
end
if nargin < 2,
    erh = warnrdlg('Input sampling rate not supplied! Defaulting to unity',...
        'SCANLINERESAMPLE : WARNING');
    uiwait(erh)
    dx = 1;
end

[ns,ntr]=size(d);
answer = inputdlg(['INPUT number of traces is ' num2str(ntr) ...
    '. Give OUTPUT Number of Traces.'],'SCANLINERESAMPLE : REQUEST',1);
if isempty(answer),              % operation canceled
    dout   = [];
    dxout  = [];
    ntrout = [];
    xout   = [];
    return
end

ntrout  = str2num(answer{1});  % output # samples per trace
r       = ntrout/ntr;          % change in sampling rate
dxout   = dx/r;                % new sampling rate
N       = 15;                  % half order of sinc summation in resample1
xout    = x0 + [0 : dxout : (ntrout-1)*dxout];

hw = helpdlg('Processing... Please wait ...');
% To comply with resample1 (working along the fast dimension), the
% radargram must be transposed
dout   = resample1(d', ntrout, ntr, N);
dout   = dout';                % restore transposed output
close(hw)

return

% END FUNCTION SCANLINERESAMPLE
