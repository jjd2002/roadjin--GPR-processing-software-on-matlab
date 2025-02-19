function x = yxtoxy(nx,dx,fx,y,ny,dy,fy,xylo,xyhi)
%
% YXTOXY : Compute a regularly-sampled, monotonically increasing function
%          x(y) from a regularly-sampled, monotonically increasing function
%          y(x) by  inverse linear interpolation.
%
% Usage  : x = yxtoxy(nx, dx, fx, y, ny, dy, fy, xylo, xyhi)
%
% Input:
% nx		number of samples of y(x)
% dx		x sampling interval; dx>0.0 is required
% fx		first x
% y         array[nx] of y(x) values; y[0] < y[1] < ... < y[nx-1] required
% ny		number of samples of x(y)
% dy		y sampling interval; dy>0.0 is required
% fy		first y
% xylo		x value assigned to x(y) when y is less than smallest y(x)
% xyhi		x value assigned to x(y) when y is greater than largest y(x)
% Output:
% x         array[ny] of x(y) values
%
% User must ensure that:
% (1) dx>0.0 && dy>0.0
% (2) y[0] < y[1] < ... < y[nx-1]
% *************************************************************************
% Transcoded from function yxtoxy.c, which is part of the CWP library.
% For function yxtoxy.c :
% Copyright (c) Colorado School of Mines, 2002.
% Author: Dave Hale, Colorado School of Mines, 06/02/89
%
% For function yxtoxy.m :
% Copyright (c) 2005, Andreas Tzanis.
% Author: Andreas Tzanis, Department of Geophysics, University of Athens,
%         atzanis@geol.uoa.gr
% *************************************************************************

nxi = nx;  dxi = dx;  fxi = fx;   nyo = ny;
dyo = dy;  fyo = fy;  fyi = y(1); jyo = 0;
yo = fyo;
% loop over output y less than smallest input y 
while jyo < nyo,
    if yo >= fyi, break, end;
    x(jyo+1) = xylo;
    jyo = jyo+1;                                                                 
    yo = yo+dyo;
end
% loop over output y between smallest and largest input y 
if jyo == nyo-1 & yo == fyi, 
    x(jyo+1) = fxi;
    jyo = jyo+1;
    yo = yo + dyo;
end
jxi1 = 0;
jxi2 = 1;
xi1 = fxi;
while jxi2 < nxi & jyo < nyo,
    yi1 = y(jxi1+1);
    yi2 = y(jxi2+1);
    if yi1 <= yo & yo <= yi2,
        x(jyo+1) = xi1+dxi*(yo-yi1)/(yi2-yi1);
        jyo = jyo+1;
        yo = yo + dyo;
    else
        jxi1 = jxi1+1;
        jxi2 = jxi2+1;
        xi1 = xi1 + dxi;
    end
end
% loop over output y greater than largest input y 
while jyo < nyo,
    x(jyo+1) = xyhi;
    jyo = jyo+1;
end
