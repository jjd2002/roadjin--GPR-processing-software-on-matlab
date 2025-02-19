function [dout, nf, lp] = predc(d,dt)  
%
%    PREDC : Predictive deconvolution of the data passed in the 2-D matrix
%            d. The program inquires the operator length, prediction
%            distance and percent prewhitening, then designs and applies a
%            deconvolution filter to each trace separately. The entire
%            trace is used in computing the auto-correlation function. 
%            Operator (filter) length must be smaller than 1/2 trace
%            length. The lag must be smaller than operator length. Length
%            and lag are given either in actual time units (ns) or in
%            number of samples.   
%
%    Usage : [dout, nf, lp] = predc(d, dt)
%
%   Inputs : 
%        d : [ns x ntr ] matrix of the GPR data (section)
%       dt : Sampling rate in nanosecs
%
%  Outputs : 
%     dout : The [ns x ntr ] deconvolved GPR section
%       nf : Length of the prediction operator used in # samples 
%       lp : Prediction distance in # samples.
%
% Requires : checkcomma.m
%            toeplitz.m - constructs Toeplitz matrix, MATLAB core 
%            pinv.m     - pseuodoinverse, MATLAB core 
%
%   Author : Andreas Tzanis,
%            Department of Geophysics,
%            University of Athens
%            atzanis@geol.uoa.gr
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

[ns, ntr]=size(d);
% Ask how deconvolution parameters will be given
ask = questdlg('Filter parameters given in # samples or in ns?',...
    'Prediction Operator parameters', ...
    'Samples','ns','Cancel','Samples');
if strcmp(ask,'Cancel'),
    dout = [];  lp=[];  nf = [];
    return
end
% Get filter parameters
% nf, Length of the filter
% lp, Prediction distance
% mu, Prewhitening in %    
answer = inputdlg({'Prediction Operator length' ... 
    'Prediction length' ...
    '% Prewhitening'},...
    'Prediction Operator parameters',1);
if isempty(answer),
    dout = [];  lp=[];  nf = [];
    return
end
ansr = checkcomma(answer); 
nf = str2num(ansr(1,:));
lp = str2num(ansr(2,:));
mu = str2num(ansr(3,:));
if strcmp(ask,'ns'),
    nf = floor(nf/dt);
    lp = floor(lp/dt);
end

% Error trapping
if nf >= ns/2, 
    erh = errordlg(['Operator length cannot be longer than ' ...
                    'half the data length!               '], ...
                'PREDC : ERROR');
    uiwait(erh);
    dout = [];  lp=[];  nf = [];
    return
end
if lp >= nf,
    erh = errordlg(['The prediction distance cannot be longer ' ...
                    'than the length of the operator!           '], ...
                'PREDC : ERROR');
    uiwait(erh);
    dout = [];  lp=[];  nf = [];
    return
end

%%% All is OK, Process ...
hw = waitbar(0,'Deconvolving, Please wait ... ');
dout = zeros(ns,ntr);
for i = 1:ntr,
    w   = d(:,i);                                    % make w column vector
    cc  = cross(w,w,nf);                   % Cross correlation (row vector)
    R   = toeplitz(cc) + ((cc(1)*mu/100)*eye(nf));        % Toeplitz Matrix
    rhs = cross(w,w,nf+lp);                      % Right hand side vector 1
    rhs = rhs(lp+1:nf+lp)';                      % Right hand side vector 2
    %f = pinv(R)*rhs;                                              % Filter 
    f = R\rhs;                                                     % Filter 
    if lp==1; 
        f = [1; -f]; 
    else
        f = [1; zeros(lp-1,1); -f];
    end
    co = conv(f,w);                         %  Actual output
    dout(:,i) = co(1:ns);    
    waitbar(i/ntr,hw);
end
close(hw);
return

function cc = cross(x,y,lc)
% Cross-correlation function of column vector x sliding through column 
% vector y through lc lags
lx = length(x);
x  = [x; zeros(lc,1)];
ly = length(y);
cc = zeros(1,lc);
for i = 1:lc
     cc(i) = sum(x(i:i+lx-1).*y(1:ly));
end
return
