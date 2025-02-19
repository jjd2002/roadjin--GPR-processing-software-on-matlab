function cmap = twocolorsc(dmin, d1, d2, dmax, ncolors, color1, color2) 
%
% TWOCOLORSC : Make a 3-color colormap of length ncolors, scaled in the
%              range [dmin,d1] and [d2,dmax], where [dmin, dmax] are
%              respectively the maximum and minimum values of the data to
%              be displayed. The colormap comprises:
%              * RBG color1 to white reel for data in range [dmin, d1], and  
%              * White to RGB color2 reel for data in range [d2, dmax]
%              * The cutoff values d1 and d2 are determined by the user
%
%      Usage : cmap = twocolorsc(dmin, d1, d2, dmax, ncolors, color1, color2) 
%
% Author : Andreas Tzanis
%          Department of Geophysics, 
%          University of Athens
%          atzanis@geol.uoa.gr
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
%

% find indices corresponding to data values d1 and d2 
dd     = (dmax - dmin)/(ncolors-1);
cindex = dmin:dd:dmax;
i1     = find(cindex <= d1);
if length(i1)==0,
    i1 = 1;
else
    i1 = i1(length(i1))+1;
end
i2     = find(cindex >= d2);
i2     = i2(1);
% make red color reel
c1      = ones(ncolors,1);
dc1a    = (1-color1(1))/(i1-1); 
if dc1a ~=0,
    c1a     = color1(1):dc1a:1;
    c1(1:i1)= c1a(1:i1)';
end
dc1b    = (1-color2(1))/((ncolors-i2));
if dc1b ~=0,
    c1b     = 1:-dc1b:color2(1);
    c1(i2:ncolors)= c1b(1:length(c1b))';
end
% make green color reel
c2      = ones(ncolors,1);
dc2a    = (1-color1(2))/(i1-1); 
if dc2a ~=0,
    c2a     = color1(2):dc2a:1;
    c2(1:i1)= c2a(1:i1)';
end
dc2b    = (1-color2(2))/((ncolors-i2));
if dc2b ~=0,
    c2b     = 1:-dc2b:color2(2);
    c2(i2:ncolors)= c2b(1:length(c2b))';
end
% make blue color reel
c3      = ones(ncolors,1);
dc3a    = (1-color1(3))/(i1-1); 
if dc3a ~=0,
    c3a     = color1(3):dc3a:1;
    c3(1:i1)= c3a(1:i1)';
end
dc3b    = (1-color2(3))/((ncolors-i2));
if dc3b ~= 0,
    c3b     = 1:-dc3b:color2(3);
    c3(i2:ncolors)= c3b(1:length(c3b))';
end
% make colormap
cmap = [c1 c2 c3];

return