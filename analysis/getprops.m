function properties = getprops(fc,props)
%
%   GETPROPS : Function used by BUILD2DMODEL, to set, or reset, the EM
%              properties of a structural element in the synthetic model.
%              The routine inquires for the resistivity, Relative
%              dielectric constant and relative permeability of the object
%              via a dialog box, converses with the user to optimize their
%              values and upon confirmation, calculates the quality factor
%              after Bano (1996, Geophys. J. Int., 279 - 288).
%
%      Usage : properties = getprops(fc,props)
%              properties = getprops(fc)
%
%     Inputs : 
%         fc : Antenna central frequency.
%      props : [1 x 4] vector holding the EM properties of an structural
%              object (for resetting). 
%
%    Outputs :
% properties : The EM properties of the object. These are:
%              a) relative dielectric constant,
%              b) quality factor, 
%              c) relative magnetic permeability,
%              d) resistivity 
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
w     = 2*pi*fc*1e6;
if nargin==1 || isempty(props),
    rho   = 1000;                 % default EM properties give v ~ 0.1m/ns
    K     = 9;
    M     = 1;
else
    rho = props(4);
    K   = props(1);
    M   = props(3);
end
    
reply = 'No';
while strcmp(reply,'Yes')~=1,
clb = cell(3,1);                                      
clb(1) = cellstr(['Give Resistivity']);
clb(2) = cellstr(['Give Relative Dielectric Constant ']);
clb(3) = cellstr(['Give Relative Magnetic Permeability']);
cldef(1) = cellstr(num2str(rho));
cldef(2) = cellstr(num2str(K));
cldef(3) = cellstr(num2str(M));
answer = inputdlg(clb,'Give body properties',1,cldef);  
if isempty(answer), 
    properties = []; 
    return; 
end;  
lb = char(answer);                                      
for i=1:3
    comma = findstr(lb(i,:),',');
    if ~isempty(comma), 
        for j=1:length(comma); 
            lb(i,comma) = '.'; 
        end;
        clear comma
    end
end
rho   = str2num(lb(1,:));                                  
sigma = 1/rho;
K     = str2num(lb(2,:));
M     = str2num(lb(3,:));
e     = K*e0; 
m     = M*m0;
k = sqrt( -w*w*e*m + sqrt(-1)*w*m/rho);                  % wavenumber
v = w/imag(k)/1e9;                                       % phase velocity
str = ['This will give a velocity of ' num2str(v) ' m/ns. OK?'];
reply = questdlg(str);
     if strcmp(reply,'Cancel')==1,
         properties = []; 
         return
     end
end
% Now that we're satisfied, generate quality factor as well
sigmastar = sigma + sqrt(-1)*w*e; 
estar     = -sqrt(-1)*sigmastar/w; 
losstan   = (sigma -w*imag(estar))/(w*real(estar)); 
Q         = 1/losstan;
properties = [K Q M rho];
return
%This is approximately air:  rho = 70000; sigma=1/rho; e=1*e0; m=1*m0;
%sigma=0.00001390758 and e=e0 will give a Q of approx. 1000;
%This gives a velocity of ~ 0.1 m/ns: %rho = 700; sigma=1/720; m=m0; e=9e0;
%Non-dispersive term of velocity given K and Q, after Bano (1996):
%v = 1/(sqrt(m*e)*cos((pi/4)*(1-((2/pi)*atan(Q)))))/1e9;