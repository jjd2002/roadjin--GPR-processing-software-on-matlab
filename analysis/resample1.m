  
function y = resample1(x,p,q,order)
%
%RESAMPLE1 : Change the sample rate of x by a factor of p/q. Note that p
%and q do % not need to be integers since this routine does not use a
%polyphase rate change algorithm, but instead uses bandlimited
%interpolation, wherein the continuous time signal is estimated by summing
%the sinc  functions of the nearest neighbouring points up to distance d.
%This is discussed in: J.O. Smith and P. Gossett (1984). A flexible
%sampling-rate conversion method. In ICASSP-84, Volume II, pp.
%19.4.1-19.4.2. New York: IEEE Press. Also see the authors page at: 
%http://www-ccrma.stanford.edu/~jos/resample/ 
%
%  Usage : y = resample1(x,p,q,order)
%
% Inputs :  x  is the 2-D data matrix to be resampled
%         p,q  determine the new sampling rate
%              order  is the half order of sinc summation 
%
%Outputs :  y  is the resampled data matrix
%
% Copyright (C) 2000, Paul Kienzle.
% Modified for MATLAB / MATGPR by: 
% Andreas Tzanis, 
% Department of Geophysics,
% University of Athens,
% atzanis@geol.uoa.gr
%
% => The program can now resample a rectangular matrix along any one of its
% dimensions 
%

if (nargin < 2 | nargin > 5)
    disp(' usage : y=resample(x,p,q,order)')
    return
end

if (nargin < 3), 
    q=1; 
end
if (nargin < 4), 
    order = 5; 
end

%  chain to decimate/interpolate if appropriate
if p==1 & q==fix(q),
    y=decimate(x,q); 
    return;
elseif q==1 & p==fix(p)
    y=interp(x,q); 
    return;
end
 
[rowsx,colsx]=size(x);
transpose = (rowsx == 1);   
if transpose, 
    x = x.'; 
end
r=p/q;
% if rate reduction, apply antialiasing filter first
if (r < 1),
    b = fir_f1(2*order+1, r);
    X = fft(x);
    b1 = [b; zeros(rowsx-length(b),1)];
    B=fft(b1) * ones(1,colsx);
    X = X.*B;
    X = X.*conj(B);
    x = real(ifft(X));
end
%% Determine the new sampling times, and their distance to the old
%% ones.  Note that the new series should be the maximum that can
%% be contained in the old series without going over the time
%% allotted to the old series.  In short, you have to go a little
%% beyond the last sample of the old series if your new sampling
%% rate is higher.
t   = [1: 1/r: (rowsx+1) - 1/r]';   % the sampling points of the new series
idx = fix(t);                       % the nearest old point
t   = t-idx;                        % distance to the nearest old point
%% generate the new series by summing the sinc functions of the
%% nearest neighbour points implicit in the continuous time
%% expansion of the old series.  This new series is truncated
%% to +/- order nearest neighbours.  For convenience, the original
%% series is zero-padded before and after, implicitly setting the
%% neighbours at the start of the signal to zero.
[rowsx,colsx]=size(x);
x = [zeros(order,colsx) ; x ; zeros(order,colsx)];
y = zeros(length(idx),colsx);              % the new series
for i=-order:order
    w = sinc(t-i).*(0.5+0.5*cos(pi*(t-i)/(order+0.5))); % hanning window
    y = y + x(idx+i+order,:).*w(:,ones(size(x,2),1));
end
if transpose, 
    y=y.'; 
end
% End Function RESAMPLE1
%
% Speech example
%  [x, fs] = auload(file_in_loadpath("sample.wav"));
%  sound(resample(x,16000,fs), 16000);  % resample at 16 kHz
% 
% Example from interp1
%  xf=0:0.05:10.95; yf = sin(2*pi*xf/5);
%  xp=0:10;         yp = sin(2*pi*xp/5);
%  r = resample1(yp,xp(2),xf(2));
%  plot(xf,yf,';original;',xf,r,';resample;',xp,yp,'*;;');
% Note that resample computes all samples up to but not including time
% n+1. If you are increasing the sample rate, this means that it will
% generate samples beyond the end of the time range of the original
% signal. That is why xf must goes all the way to 10.95 in the example.
%
