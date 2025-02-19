function slice3d(D3D)
%
% SLICE3D : Generates 3-D Slice Plots of GPR Profile Data With or Without
%           Translucence
%
%Comments : Presupposes the following Co-ordinate system:
%            ^
%            |
%            |...........................................>
%            |    Direction of survey (scan) lines 
%     Y-axis |...........................................>
%            |
%            |...........................................>
%            |
%            |____________________________________________>
%                      Scan Axis (X-axis)
%
%    ==> The survey is conducted along parallel lines (profiles). Each
%        3-D slice corresponds to a GPR profile
%    ==> The X-axis is the longitudinal direction of the profiles 
%    ==> The Y-axis os the meridional direction, i.e. perpendicular to the
%        profiles.
%    ==> Input data structure D3D generated with CONCAT3D_2 and contains:
%        D3D.x : 3-D matrix with the X-coordinates of GPR traces in the
%                profiles  
%        D3D.y : 3-D matrix with the Y-coordinates of GPR traces in the
%                profiles  
%        D3D.z : 3-D matrix with the Z-coordinates of data samples. This
%                dimension can be either time or depth
%        D3D.d : The 3-D GPR data volume with coordinates [D3D.x, D3D.y,
%                D3D.z]
%      3D.zlab : The vertical axis label, either time or depth
% 
% Usage : slice3d(D3D)
%
%Author : Andreas Tzanis,
%         Department of Geophysics, 
%         University of Athens
%         atzanis@geol.uoa.gr
%
% Copyright (C) 2008, Andreas Tzanis. All rights reserved.
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

% Initialize 
figure('name','3D Slice Display','tag','slice3dfigure',...
    'menubar','none', 'numbertitle','off', ...
    'position',[80   300   700   500], ...
    'DeleteFcn',['matgprwindows(''updateremove'', ' ...
                 'findobj(''tag'',''slice3dfigure'')); ' ...
                 'delete(findobj(''tag'',''slice3dcontrols'')); ' ...
                 'clear hs zs vpos xzax yax asp;']); 
% Setup he MATGPR figure manipuation utilities
figuretools; 
imagecolors;
matgprwindows('setup'); 
matgprwindows('updateadd');

% Basic plot 
x   = D3D.x;
y   = D3D.y;
z   = D3D.z;
d3d = D3D.d;
% Fix axes limits to data limits
xlim = [min(x(1,:,1)) max(x(1,:,1))];
ylim = [min(y(:,1,1)) max(y(:,1,1))];
zlim = [min(z(1,1,:)) max(z(1,1,:))];
hs = slice(x,y,z,d3d,[],y(:,1,1),[]);
% Make slice plot
set(hs,'FaceColor','flat','EdgeColor','none','tag','h3dslices');
set(gca,'tag','3dslicesaxes', ...
    'zdir','reverse', ...
    'ydir','reverse', ...
    'nextplot', 'add', ...
    'xlim', xlim, ...
    'ylim', ylim, ...
    'zlim', zlim );
box on;                            % Puts box around data for visualization
grid off;
xlabel('Scan Axis in m (x-axis)');
ylabel('Y-Axis (m)');
zlabel(D3D.zlab);
% Set axis aspect ratio
daspect([1 1 10])

% Initialize gui controls
slice3d_gui(x,y,z,d3d)

return
