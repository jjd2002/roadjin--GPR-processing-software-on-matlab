function  dwow = dewow( d )
%
%    DEWOW  : Eliminate wow by applying a zero-phase high pass FIR filter
%    frequencies with cutoff frequency at 2% of the Nyquist  
%
%     Usage : dwow = dewow( d )
%
%    Inputs : 
%        d  : 2-D array of distance vs time input data;  
%   Outputs : 
%      dwow : 2-D array of filtered data
%
% Requires  : fir_f1.m, fir_f2.m 
%
%   Author  : Andreas Tzanis,
%             Department of Geophysics, 
%             University of Athens
%             atzanis@geiol.uoa.gr
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
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

%%% ns and ntr is Samples_per_scan and No_traces respectively
[ns, ntr] = size(d);
%%% Set cutoff at 2% of the Nyquist
Ws = 0.02;
%%% Make filter long for stability 
nf = floor(0.9*ns);
if mod(nf,2) ~= 0,
    nf = nf+1;
end
b  = fir_f1(nf, Ws, 'high');
b = [b; zeros(ns-length(b),1); zeros(ns,1)];
%%%% Filter in the frequency domain and twice to preserve phase
hw = helpdlg('DEWOW: Processing... Please wait ...');
H = fft(b); 
H     = H * ones(1, ntr);         % cast filter into a matrix operator
fttr  = fft([d; zeros(ns,ntr)]);  % transform data to frequecy domain
fttr  = fttr.*H;                  % filter forwards
fttr  = fttr.*conj(H);            % filter backwards
dwow  = real(ifft(fttr));         % invert to time domain
dwow  = dwow(1:ns,:);

close(hw);
return
