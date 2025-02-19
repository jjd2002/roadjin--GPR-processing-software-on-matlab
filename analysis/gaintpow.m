function [dtpow, scale, pow] = gaintpow( d, traveltime );
%
%  GAINTPOW : * Multiplies the columns (traces) of the input data matrix "d"
%               with a gain function of the form g(t) = scale * t^pow. 
%             * The exponent "pow" is requested from the user. However, the
%               program uses function "getattenuation.m" to compute the
%               attenuation characteristics of the data and determine a best
%               fitting model, suggesting the resulting "pow" to the
%               user. This is, usually, a nearly optimal solution, but 
%               experimentation may provide a result more suitable for some
%               particular needs / requirements.
%             * The attenuation characteristics and the best fitting
%               power-law model can be visualized by running
%               "getattenuation" as:  getattenuation(d, traveltime);
%             * The "scale" is computed automatically, so as to preserve the
%               amplitude range values of the input data (compensate for the
%               effects of "t^pow").
%
%     Usage : [dtpow, scale, pow]  = gaintpow( d, dt );
%   
%    Inputs :
%         d : [ns x ntr) input data array (GPR section).
%traveltime : [1 x ns] two-way traveltime vector.
%
%   Outputs : 
%     dtpow : The output (amplified) data array
%     scale : The scale factor in the equation g(t) = scale * t^pow.
%       pow : The exponent in the equation g(t) = scale * t^pow.
%
% Requires : getattenuation.m
%
%   Author : Andreas Tzanis
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


[ns,ntr]   = size(d);
p = getattenuation(d, traveltime, 'no_plot');
% Set defaults
pow   = abs(p(2))-1;
default = {num2str(pow)};
prompt = [ ...
'请提供指数：该值建议低于由数据计算出的最佳拟合幂律衰减模型的结果。' ];
ask = inputdlg(prompt, 'Gain function g(t) = scale*t^pow', 1, default);
if isempty(ask), 
    dtpow = [];  scale=[];  pow=[];
    return; 
end;
lb = char(ask);                   % check for commas and replace with dots
comma = findstr(lb(1,:),',');
if ~isempty(comma),
    for j=1:length(comma); 
         lb(i,comma) = '.'; 
    end;
end
pow   = str2num(lb(1,:));

g     = (traveltime.^pow)' * ones(1,ntr);
dtpow = d.*g;

% scale factor
s2 = (max(max(dtpow))-min(min(dtpow)));
s1 = (max(max(d))-min(min(d)));
scale = s1/s2;
dtpow = scale*dtpow;
return
