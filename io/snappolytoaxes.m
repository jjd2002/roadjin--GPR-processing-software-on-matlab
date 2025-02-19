function [x, z] = snappolytoaxes(xi, zi, xl, zl);
%
% SNAPPOLYTOAXES : Clip a polygonal area to axes limits. For a given
%                  polygon with vertices [x,z] and an axes object with
%                  limits xl = [xlow, xhigh] and zl = [zlow, zhigh], the
%                  routine checks for the indices ii and jj of outlier
%                  vertices: x(ii)< xl < x(jj) & z(ii)< zl < z(jj). 
%                  Those found coming on or going off are snapped onto the
%                  axes limits, by linear interpolation. The rest are
%                  snapped by brute force.  
%              ==> This program is part of the MATGPR modelling suite.
%              ==> Will never claim credit for intelligent programming
%                  solutions!
%
%          Usage : [x, z] = snappolytoaxes(x, z, xl, zl)
%
%         Inputs : 
%         xi, zi : Column vectors containing the polygon vertices
%         xl, zl : The horizontal and vetical limits of the axes object
%                  containing the polygon.
%        Outputs : 
%           x, z : Column vectors containing the snapped polygon vertices
%
%         Author : Andreas Tzanis
%                  Department of Geophysics, 
%                  University of Athens
%                  atzanis@geol.uoa.gr
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

% Close polygon for snapping
x = [xi; xi(1)];
z = [zi; zi(1)];
np = length(x);
% begin
for i = 2:np,
    if x(i) > xl(1) && x(i-1) < xl(1), 
    slope = (z(i) - z(i-1))/(x(i) - x(i-1));
        z(i-1) = z(i) + slope*(xl(1)-x(i));
        x(i-1) = xl(1);
    end
    if x(i-1) > xl(1) && x(i) < xl(1),
    slope = (z(i-1) - z(i))/(x(i-1) - x(i));
        z(i) = z(i-1) + slope*(xl(1) - x(i-1));
        x(i) = xl(1);
    end
end
ii    = find(x < xl(1));
x(ii) = xl(1);
% check right
for i = 2:np,
    if x(i) > xl(2) && x(i-1) < xl(2),
        slope = (z(i) - z(i-1))/(x(i) - x(i-1));
        z(i) = z(i) + slope*(xl(2)-x(i));
        x(i) = xl(2);
    end
    if x(i) < xl(2) && x(i-1) > xl(2),
        slope = (z(i-1) - z(i))/(x(i-1) - x(i));
        z(i-1) = z(i-1) + slope*(xl(2)-x(i-1));
        x(i-1) = xl(2);
    end
end
ii    = find(x > xl(2));
x(ii) = xl(2);
% check top
for i = 2:np,
    if z(i) < zl(1) && z(i-1) > zl(1),
        den = (x(i) - x(i-1));
        if den ~=0,
            slope = (z(i) - z(i-1))/den;
            z0 = z(i) + slope*(xl(1)-x(i));
            x(i) = (zl(1) - z0)/slope;
            z(i) = zl(1);
        else
            z(i) = zl(1);
        end
    end
    if z(i) > zl(1) && z(i-1) < zl(1),
        den = (x(i-1) - x(i));
        if den~=0,
            slope = (z(i-1) - z(i))/den;
            z0 = z(i-1) + slope*(xl(1)-x(i-1));
            x(i-1) = (zl(1) - z0)/slope;
            z(i-1) = zl(1);
        else
            z(i-1) = zl(1);
        end
    end
end
ii = find(z < zl(1));
z(ii) = zl(1);
% check bottom
for i = 2:np,
    if z(i) < zl(2) && z(i-1) > zl(2),
        den = (x(i) - x(i-1));
        if den ~=0, 
            slope = (z(i) - z(i-1))/den;
            z0 = z(i) + slope*(xl(1)-x(i));
            x(i-1) = (zl(2) - z0)/slope;
            z(i-1) = zl(2);
        else
            z(i-1) = zl(2);
        end
    end
    if z(i) > zl(2) && z(i-1) < zl(2),
        den = (x(i-1) - x(i));
        if den ~= 0, 
            slope = (z(i-1) - z(i))/den;
            z0 = z(i-1) + slope*(xl(1)-x(i-1));
            x(i) = (zl(2) - z0)/slope;
            z(i) = zl(2);
        else
            z(i) = zl(2);
        end
    end
end
ii = find(z > zl(2));
z(ii) = zl(2);
%%% Finished - reset polygon to oroginal size
x = x(1:np-1);
z = z(1:np-1);
return