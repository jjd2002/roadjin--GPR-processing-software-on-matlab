function [dout, dtout, nsout] = timeresample(d, dt)

% TIMERESAMPLE : Driver for routine RESAMPLE1.M - resamples the traces of
%                the 2-D GPR data matrix d at a higher or lower rate. 
%
%        Usage : [dout, dtout, nsout] = timeresample(d, dt)
%
%       INPUTS : 
%       d      : input 2-D radargram
%       dt     : input sampling rate
%
%      OUTPUTS : 
%       dout   : the resampled radargram
%       dtout  : the sampling rate in dout.
%       nsout  : the # samples in dout
%
%     REQUIRES : RESAMPLE1.M
%
%           Author : Andreas Tzanis, 
%                    Department of Geophysics, 
%                    University of Athens
%                    atzanis@geol.uoa.gr
%                   (C) 2005, Andreas Tzanis, all rights reserved
%
if nargin < 2,
    erh = warnrdlg(['Input sampling rate not supplied! ' ...
        'Defaulting to unity'],'TIMERESAMPLE : WARNING');
    uiwait(erh)
    dt = 1;
end
[rowsd,colsd]=size(d);
answer = inputdlg(['INPUT Samples per scan are ' num2str(rowsd) ...
    '. Give OUTPUT Samples per scan.'], 'TIMERESAMPLE : REQUEST',1);
if isempty(answer),                 % operation canceled
    dout  = [];
    dtout = [];
    nsout = [];
    return
end
nsout  = str2num(answer{1});   % output # samples per trace
r      = nsout/rowsd;          % change in sampling rate
dtout  = dt/r;                 % new sampling rate
N      = 15;                   % half order of sinc summation in resample1

dout   = resample1(d, nsout, rowsd, N);

return
% END FUNCTION TIMERESAMPLE
