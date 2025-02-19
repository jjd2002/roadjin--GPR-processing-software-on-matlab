function ac = removedc( d );
%
%  REMOVEDC : Subtracts DC component (mean) from each individual trace
%  (column) of the 2-D GPR data matrix. 
%
%   Usage : ac = removedc( d );
%
%   Input :  d,  the 2-D GPR section
%
%  Output : ac,  the reduced GPR section
%
%  Author : Andreas Tzanis,
%           Department of Geophysics, 
%           University of Athens
%           atzanis@geol.uoa.gr
%
%  (C) 2006, Andreas Tzanis, all rights reserved
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

[ns,ntr] = size(d);
dc = ones(ns,1)*mean(d,1);
ac = d - dc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The old fashioned way will not do any more ...
% ac       = zeros(ns,ntr); 
% h = waitbar(0,'Removing DC ...');
% for i=1:ntr; 
%     ac(:,i) = d(:,i) - mean(d(:,i)); 
%     waitbar(i/ntr,h);
% end;
% close(h);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return