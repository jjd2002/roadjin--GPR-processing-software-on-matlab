function [scanline,tt2w] = hyperbola(xstart,xstop,TxRx,x,d,r,v)   
%
% HYPERBOLA : Calculate a difraction front hyperbola for a simple object,
%             - point diffractor, (quasi)cylinder or )quasi)sphere -
%             assuming that it is embedded in a non-dispersive uniform
%             halfspace (constant velocity).
%
%     Usage : [scanline,tt2w] = hyperbola(xstart,xstop,TxRx,x,d,r,v)   
%
%   Inputs  : 
%     xstart: First point on the scan line to (in meters);
%     xstop : Last point of the scan line (meters);
%      TxRx : Source - Receiver offset; 
%         x : Location of the centre of the object relative to
%             xstart (in m) 
%         d : Depth to the centre of the object (in m)
%         r : Object radius (in m)
%         v : Velocity (in m/ns)
%
%   Outputs : 
%  scanline : Coordinates of the midpoint between Tx and Rx offset 
%      tt2w : 2-way traveltime (hyperbola)
%
%  Author  : A. Tzanis
%            Department of Geophysics,
%            University of Athens.
%            atzanis@geol.uoa.gr
%
%  (C) 2005, Andreas Tzanis, all rights reserved
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

TxRx_midpoint = TxRx/2.0;
xstep = (xstop - xstart) / 200.0;
scanline = xstart : xstep : xstop;
Tx_location = scanline - TxRx_midpoint;
Rx_location = scanline + TxRx_midpoint;
Tx2target = sqrt((x - Tx_location).*(x - Tx_location) + d*d);
Target2Rx = sqrt((x - Rx_location).*(x - Rx_location) + d*d);
Tx2target = Tx2target - r;
Target2Rx = Target2Rx - r;
Tx2Rx = Tx2target + Target2Rx;
tt2w = Tx2Rx / v;
return

