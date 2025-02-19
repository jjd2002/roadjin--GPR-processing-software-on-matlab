function isosurf(D3D)
%
% ISOSURF : Generates isometric surfaces of equal GPR signal intensity 
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
%    ==> The survey is conducted along parallel scan lines (profiles)
%    ==> The X-axis is the longitudinal direction of the profiles 
%    ==> The Y-axis os the meridional direction, i.e. perpendicular to the
%        profiles.
%    ==> Input data structure D3D generated with CONCAT3D_2 and contains:
%        D3D.x : 3-D matrix with the X-coordinates of GPR traces along the
%                profiles  
%        D3D.y : 3-D matrix with the Y-coordinates of GPR traces across the
%                profiles  
%        D3D.z : 3-D matrix with the Z-coordinates of data samples. This
%                dimension can be either time or depth
%        D3D.d : The 3-D GPR data volume with coordinates [D3D.x, D3D.y,
%                D3D.z]
%      3D.zlab : The vertical axis label, either time or depth
% 
% Usage : isosurf(D3D)
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

figure('name','Isosurface Display','tag','isofigure',...
    'menubar','none', 'numbertitle','off', ...
    'position',[80   300   700   500], ...
    'DeleteFcn',['matgprwindows(''updateremove'', ' ...
                 'findobj(''tag'',''isofigure'')); ' ...
                 'delete(findobj(''tag'',''iso3dcontrols'')); '...
                 'clear rdf rds rdv lh; ']); 
% Setup he MATGPR figure manipuation utilities
figuretools;
matgprwindows('setup'); 
matgprwindows('updateadd');

% Basic Plot
x   = D3D.x;
y   = D3D.y;
z   = D3D.z;
d3d = D3D.d;
isovalue  = 0.55*max(max(max(d3d)));
% Create Isosurface
[rdf,rdv] = isosurface(x,y,z,d3d,isovalue);
rds = patch('Faces',rdf,'Vertices',rdv);
set(rds,'tag','isosurface','edgecolor','none','facecolor',[0.58,0.39,0.39]);
% Surface caps to plug holes
%[cf,cv] = isocaps(x,y,z,d3d,isovalue,'above');
%rdc = patch('Faces',cf,'Vertices',cv);
%set(rdc,'tag','isocaps','edgecolor','none','facecolor',[0.58,0.39,0.39])
% Surface normals for smoother presentation
isonormals(x,y,z,d3d,rdf);

% Fix axes limits to data limits
xlim = [min(x(1,:,1)) max(x(1,:,1))];
ylim = [min(y(:,1,1)) max(y(:,1,1))];
zlim = [min(z(1,1,:)) max(z(1,1,:))];
set(gca,'tag','isoaxis',...
    'zdir','reverse',...
    'ydir','reverse', ...
    'xlim', xlim, ...
    'ylim', ylim, ...
    'zlim', zlim );
view(3);
box on;
grid on;
xlabel('Scan Axis in m (x-axis)');
ylabel('Y-Axis (m)');
zlabel(D3D.zlab);
% Set axis aspect ratio
daspect([1 1 10])

% Create light source
lh = camlight('headlight');
lighting gouraud
set(lh,'tag','isolight');

% Initialize GUI controls
isosurf_gui(x,y,z,d3d);

return