function b = fir_f1(n, w, ftype, window, scale)
%
% Produce an order n FIR filter with the given frequency cutoff,
% returning the n+1 filter coefficients in b.  
%
% usage: b = fir_f1(n, w [, type] [, window] [, noscale])
%
% n: order of the filter (1 less than the length of the filter)
% w: band edges
%    strictly increasing vector in range [0, 1]
%    singleton for highpass or lowpass, vector pair for bandpass or
%    bandstop, or vector for alternating pass/stop filter.
% type: choose between pass and stop bands
%    'high' for highpass filter, cutoff at w
%    'stop' for bandstop filter, edges at w = [lo, hi]
%    'DC-0' for bandstop as first band of multiband filter
%    'DC-1' for bandpass as first band of multiband filter
% window: smoothing window
%    defaults to hamming(n+1) row vector
%    returned filter is the same shape as the smoothing window
% noscale: choose whether to normalize or not
%    'scale': set the magnitude of the center of the first passband to 1
%    'noscale': don't normalize
%
% To apply the filter, use the return vector b:
%       y=filter(b,1,x);
%
% Examples:
%   freqz(fir_f1(40,0.3));
%   freqz(fir_f1(15,[0.2, 0.5], 'stop'));  % note the zero-crossing at 0.1
%   freqz(fir_f1(15,[0.2, 0.5], 'stop', 'noscale'));
%
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
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
%
% => Original function name was fir1.m, built for OCTAVE. 
% => Adapted for MATLAB/MATGPR by Andreas Tzanis, April 2004.
% => Renamed to fir_f1.m so as not to conflict with MATLAB’s own fir1
%    function (Signal Analysis toolbox).
%

if nargin < 2 || nargin > 5,
    erh = helpdlg('Usage :b = fir_f1(n, w [, type] [, window] [, noscale])',...
        'FIR_F1 : HELP');
    uiwait(erh)
    return
end

% interpret arguments
if nargin==2,
    ftype=''; window=[]; scale=[];
elseif nargin==3,
    window=[]; scale=[];
    if ~ischar(ftype), 
        window=ftype; 
        ftype=''; 
    end
elseif nargin==4,
    scale=[];
    if ischar(window), 
        scale=window; 
        window=[]; 
    end
    if ~ischar(ftype), 
        window=ftype; 
        ftype=''; 
    end
end

% If single band edge, the first band defaults to a pass band
% to create a lowpass filter.  If multiple band edges, assume
% the first band is a stop band, so that the two band case defaults
% to a band pass filter.  Ick.
ftype = lower(ftype);
if isempty(ftype), 
    ftype = length(w)==1;
elseif strcmp(ftype, 'low'), 
    ftype = 1;
elseif strcmp(ftype, 'high'), 
    ftype = 0;
elseif strcmp(ftype, 'pass'), 
    ftype = 0;
elseif strcmp(ftype, 'stop'), 
    ftype = 1;
elseif strcmp(ftype, 'dc-0'), 
    ftype = 0;
elseif strcmp(ftype, 'dc-1'), 
    ftype = 1;
elseif ischar(ftype),
    erh = errordlg(['Invalid filter type ', ftype],'fir_f1 : ERROR');
    uiwait(erh)
    b = [];
    return
else
    erh = errordlg('Filter type should be a string','fir_f1 : ERROR');
    uiwait(erh)
    b = [];
    return
end

% scale the magnitude by default
if isempty(scale) || strcmp(scale, 'scale'),
    scale = 1;
elseif strcmp(scale, 'noscale'),
    scale=0;
else
    erh = errordlg('Scale must be ''scale'' or ''noscale''',...
        'fir_f1 : ERROR');
    uiwait(erh)
    b = [];
    return;
end

% use fir_f2 default filter
if isempty(window) && isempty(scale), 
    window = []; 
end

% build response function according to fir_f2 requirements
bands = length(w)+1;
f = zeros(1,2*bands);
f(1) = 0; f(2*bands)=1;
f(2:2:2*bands-1) = w;
f(3:2:2*bands-1) = w;
m = zeros(1,2*bands);
m(1:2:2*bands) = rem([1:bands]-(1-ftype),2);
m(2:2:2*bands) = m(1:2:2*bands);

% Increment the order if the final band is a pass band.  Something
% about having a nyquist frequency of zero causing problems.
if rem(n,2)==1 && m(2*bands)==1, 
    disp(['fir_f1 : Filter order (n) must be even for highpass ' ... 
        'and bandstop filters. Incrementing.']);
    n=n+1; 
    if ~isempty(window), 
        rowsw = size(window,1);
        if rowsw == 1,
            window = [window(1:n/2), window(n/2:n-1)];
        else
            window = [window(1:n/2); window(n/2:n-1)];
        end
    end
end

% compute the filter
b = fir_f2(n, f, m, 2*2^nextpow2(n), 2, window);
if isempty(b),
    return
end
% normalize filter magnitude
if scale == 1,
    if m(1) == 1,                 % find the middle of the first band edge
        w_o = (f(2)-f(1))/2;
    else 
        w_o = f(3) + (f(4)-f(3))/2;
    end
    renorm = 1/abs(polyval(b, exp(-1i*pi*w_o)));     % compute |h(w_o)|^-1
    b = renorm*b;                                    % normalize the filter
end
%END FUNCTION fir_f1
%
% Examples:
%   freqz(fir_f1(40,0.3));
%   freqz(fir_f1(15,[0.2, 0.5], 'stop'));  % note the zero-crossing at 0.1
%   freqz(fir_f1(15,[0.2, 0.5], 'stop', 'noscale'));
%