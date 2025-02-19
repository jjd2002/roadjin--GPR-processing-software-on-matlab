function [v, Q] = vcalc()
%
% VCALC : Velocity Calculator. Enquires the nominal frequency of the GPR
% antenna, the resistivity, relative dielectric constant and relative
% permeability of the medium and calculates the non-dispersive term of the
% phase velocity and the quality factor after Bano (1996, Geophys. J. Int.,
% 279 - 288).  
%
%   Usage : [v, Q ] = vcalc;
%
% Outputs :
%       v : The non-dispersive term of the phase velocity for the given
%           medium and frequency.
%       Q : The quality factor for the given medium and frequency.
%
%  Author : Andreas Tzanis
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

e0    = 8.854e-12;
m0    = 4*pi*1e-7;
fc    = 250;
rho   = 1000;                 % default EM properties give v ~ 0.1m/ns
K     = 9;
M     = 1;
    
reply = 'No';
while strcmp(reply,'OK')~=1,
    clb = cell(4,1);                                      
    clb(1) = cellstr('Give Antenna Central Frequency in MHz');
    clb(2) = cellstr('Give Resistivity');
    clb(3) = cellstr('Give Relative Dielectric Constant ');
    clb(4) = cellstr('Give Relative Magnetic Permeability');
    cldef(1) = cellstr(num2str(fc));
    cldef(2) = cellstr(num2str(rho));
    cldef(3) = cellstr(num2str(K));
    cldef(4) = cellstr(num2str(M));
    answer = inputdlg(clb,'VCALC: Give EM properties',1,cldef);  
    if isempty(answer), 
        v=[]; Q=[];  
        return; 
    end;  
    lb = char(answer);                                      
    for i=1:4
        comma = findstr(lb(i,:),',');
        if ~isempty(comma), 
            for j=1:length(comma); 
                lb(i,comma) = '.'; 
            end;
            clear comma
        end
    end
    fc    = str2num(lb(1,:));
    w     = 2*pi*fc*1e6;
    rho   = str2num(lb(2,:));                                  
    sigma = 1/rho;
    K     = str2num(lb(3,:));
    M     = str2num(lb(4,:));
    e     = K*e0; 
    m     = M*m0;
    k = sqrt( -w*w*e*m + sqrt(-1)*w*m/rho);               % wavenumber
    v = w/imag(k)/1e9;                                    % phase velocity
    % Calculate quality factor as well
    sigmastar = sigma + sqrt(-1)*w*e; 
    estar     = -sqrt(-1)*sigmastar/w; 
    losstan   = (sigma -w*imag(estar))/(w*real(estar)); 
    Q         = 1/losstan;
    str1 = ['This will give a velocity of ' num2str(v) ' m/ns'];
    str2 = ['and a Quality factor of ' num2str(Q)];
    reply = questdlg(str2mat(str1,str2),'ANSWER','OK','More','OK');
    if strcmp(reply,'Cancel')==1,
        v = []; 
        Q = [];
        return
    end
end
return
%This is approximately air:  rho = 70000; sigma=1/rho; e=1*e0; m=1*m0;
%sigma=0.00001390758 and e=e0 will give a Q of approx. 1000;
%This gives a velocity of ~ 0.1 m/ns: %rho = 700; sigma=1/720; m=m0; e=9e0;
%Non-dispersive term of velocity given K and Q, after Bano (1996):
%v = 1/(sqrt(m*e)*cos((pi/4)*(1-((2/pi)*atan(Q)))))/1e9;
