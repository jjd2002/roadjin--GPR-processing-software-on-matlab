function b = fir_f2(n, f, m, grid_n, ramp_n, window)
%
% Produce an FIR filter of order n with arbitrary frequency response, 
% returning the n+1 filter coefficients in b.  
%
% usage: b = fir_f2(n, f, m [, grid_n [, ramp_n]] [, window])
%
% n: order of the filter (1 less than the length of the filter)
% f: frequency at band edges
%    f is a vector of nondecreasing elements in [0,1]
%    the first element must be 0 and the last element must be 1
%    if elements are identical, it indicates a jump in freq. response
% m: magnitude at band edges
%    m is a vector of length(f)
% grid_n: length of ideal frequency response function
%    defaults to 512, should be a power of 2 bigger than n
% ramp_n: transition width for jumps in filter response
%    defaults to grid_n/20; a wider ramp gives wider transitions
%    but has better stopband characteristics.
% window: smoothing window
%    defaults to hamming(n+1) row vector
%    returned filter is the same shape as the smoothing window
%
% To apply the filter, use the return vector b:
%       y=filter(b,1,x);
% Note that plot(f,m) shows target response.
%
% Example:
%   f=[0, 0.3, 0.3, 0.6, 0.6, 1]; m=[0, 0, 1, 1/2, 0, 0];
%   [h, w] = freqz(fir2(100,f,m));
%   plot(f,m,';target response;',w/pi,abs(h),';filter response;');
%
% Feb 27, 2000 PAK
%     use ramping on any transition less than ramp_n units
%     use 2^nextpow2(n+1) for expanded grid size if grid is too small
% 2001-01-30 PAK
%     set default ramp length to grid_n/20 (i.e., pi/20 radians)
%     use interp1 to interpolate the grid points
%     better(?) handling of 0 and pi frequency points.
%     added some demos
%
% Copyright (C) 2000 Paul Kienzle
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA
%
% => Original function name was fir2.m, , built forOCTAVE. 
% => Adapted for MATLAB/MATGPR by Andreas Tzanis, April 2004.
% => Renamed to fir_f2.m so as not to conflict with MATLAB’s own fir2 
%    function (Signal Analysis toolbox).
%
if nargin < 3 || nargin > 6
    erh = helpdlg(['Usage : b = fir_f2(n, f, m [, grid_n ' ...
        '[, ramp_n]] [, window])'], 'FIR_F2 : HELP');
    uiwait(erh);
    b = [];
    return
end
% verify frequency and magnitude vectors are reasonable
t = length(f);
if t<2 || f(1)~=0 || f(t)~=1 || any(diff(f)<0),
    erh = errordlg(['Frequency must be nondecreasing starting from 0' ...
        ' and ending at 1'],'FIR_F2 : ERROR');
    uiwait(erh);
    b = [];
    return
end
if t ~= length(m),
    erh = errordlg(['Frequency and magnitude vectors must be ' ...
        'of the same length'],'FIR_F2 : ERROR');
    uiwait(erh);
    b = [];
    return
end
% find the grid spacing and ramp width
if (nargin>4 && length(grid_n)>1) || (nargin>5 ...
        && (length(grid_n)>1 || length(ramp_n)>1)),
    erh = errordlg(['Grid size (grid_n) and transition band size ' ...
        '(ramp_n) must be integers'],'FIR_F2 : ERROR');
    uiwait(erh);
    b = [];
    return
end
if nargin < 4, 
    grid_n=512; 
end
if nargin < 5, 
    ramp_n=grid_n/20; 
end
% find the window parameter, or default to hamming
w=[];
if length(grid_n)>1, 
    w=grid_n; 
    grid_n=512; 
end
if length(ramp_n)>1, 
    w=ramp_n; 
    ramp_n=grid_n/20; 
end
if nargin < 6, 
    window=w; 
end
if isempty(window), 
    if exist('hamming.m','file') == 2,
        window=hamming(n+1); 
    else
        window=hamming1(n+1); 
    end
end
if length(window) ~= n+1, 
    erh = errordlg('Window size must be of length n+1','FIR_F2 : ERROR'); 
    uiwait(erh);
    b = [];
    return
end

% make sure grid is big enough for the window
if 2*grid_n < n+1, 
    grid_n = 2^nextpow2(n+1); 
end
% Apply ramps to discontinuities
if (ramp_n > 0)
    % remember original frequency points prior to applying ramps
    basef = f;
    %basem = m;

    % separate identical frequencies
    idx = find (diff(f) == 0);
    f(idx) = f(idx) - ramp_n/grid_n/2;
    f(idx+1) = f(idx+1) + ramp_n/grid_n/2;

    %  make sure the grid points stay monotonic
    idx = find (diff(f) < 0);
    f(idx)   = (basef(idx) + basef(idx+1))/2;
    f(idx+1) = (basef(idx) + basef(idx+1))/2;

    % preserve window shape even though f may have changed
    % m = interp1(basef, basem, f);

end

% interpolate between grid points
grid = interp1(f,m,linspace(0,1,grid_n+1)');

% Transform frequency response into time response and
% center the response about n/2, truncating the excess
b = ifft([grid ; grid(grid_n:-1:2)]);
mid = (n+1)/2;
b = real ([ b([2*grid_n-floor(mid)+1:2*grid_n]) ; b(1:ceil(mid)) ]);

% Multiplication in the time domain is convolution in frequency,
% so multiply by our window now to smooth the frequency response.
rowsw = size(window,1);
if rowsw > 1
    b = b .* window;
else
    b = b' .* window;
end

% END FUNCTION FIR_F2
%
%demo
% f=[0, 0.3, 0.3, 0.6, 0.6, 1]; m=[0, 0, 1, 1/2, 0, 0];
% plot(f,20*log10(m+1e-5),';target response;');
% hold on;
% [h, w] = freqz(fir_f2(50,f,m,512,0));
% plot(w/pi,20*log10(abs(h)),';filter response (ramp=0);');
% [h, w] = freqz(fir_f2(50,f,m,512,25.6));
% plot(w/pi,20*log10(abs(h)),';filter response (ramp=pi/20 rad);');
% [h, w] = freqz(fir_f2(50,f,m,512,51.2));
% plot(w/pi,20*log10(abs(h)),';filter response (ramp=pi/10 rad);');
% hold off;
%
%demo
% % Classical Jakes spectrum
% % X represents the normalized frequency from 0
% % to the maximum Doppler frequency
% asymptote = 2/3;
% X = linspace(0,asymptote-0.0001,200);
% Y = (1 - (X./asymptote).^2).^(-1/4);
%
% % The target frequency response is 0 after the asymptote
% X = [X, asymptote, 1];
% Y = [Y, 0, 0];
%
% title('Theoretical/Synthesized CLASS spectrum');
% xlabel('Normalized frequency (Fs=2)');
% ylabel('Magnitude');
%
% plot(X,Y,'b;Target spectrum;'); 
% hold on;
% [H,F]=freqz(fir_f2(20, X, Y));  
% plot(F/pi,abs(H),'c;Synthesized spectrum (n=20);');
% [H,F]=freqz(fir_f2(50, X, Y));  
% plot(F/pi,abs(H),'r;Synthesized spectrum (n=50);');
% [H,F]=freqz(fir_f2(200, X, Y)); 
% plot(F/pi,abs(H),'g;Synthesized spectrum (n=200);');
% hold off;
% xlabel(''); ylabel(''); title('');
%
