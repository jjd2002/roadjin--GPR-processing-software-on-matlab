function D3D = readD3D()
%
% READD3D : Import a 3-D data volume saved with MATGPR's generic M3D binary
%           format. For details about the M3D format refer to function
%           "saveD3D.m" 
%
% Usage   : D3D = readD3D();
%
% Output  : D3D, the 3-D data volume
%
% Author : Andreas Tzanis
%          Department of Geophysics, 
%          University of Athens
%          atzanis@geol.uoa.gr
%
% Copyright (C) 2008 Andreas Tzanis. All rights reserved.
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

global ENVAR
% Get file name
[xfname, xpname]= uigetfile( ...
    {'*.m3d','MATGPR M3D-Files (*.m3d)' }, ...
     'Import 3-D Data Volume ...',ENVAR.currentworkdir); 
if xfname == 0,                               % Cancelled
    D3D = [];
    return
end

[fid, result] = fopen([xpname xfname],'r',ENVAR.endian);
if fid == -1,
    erh = errordlg(result, 'MATGPR: ERROR'); 
    uiwait(erh);  
    D3D = [];
    return; 
end; 
nx = fread(fid,1,'int16');                   % dim of X-axis 
x1 = fread(fid,1,'float');                   % location of first trace
dx = fread(fid,1,'float');                   % trace spacing
xaxis = x1 : dx : (nx-1)*dx;
ny = fread(fid,1,'int16');                   % dim of Y-axis 
yaxis = fread(fid,ny,'float');
nz = fread(fid,1,'int16');                   % dim of Z-axis 
z1 = fread(fid,1,'float');                   % minimum Z-value
dz = fread(fid,1,'float');                   % Z- spacing
zaxis = z1 : dz : (nz-1)*dz;
d3d = zeros(ny, nx, nz);                     % 3-D Data volume
for i=1:ny,
    d3d(i,1:nx,1:nz) = fread(fid,[nx nz],'float');
end
ls = fread(fid,1,'int16');                   % dim of Z-axis label 
lab = fread(fid,ls,'uchar');                 % Z-axis label
zlab = char(lab');
fclose(fid);
[x,y,z] = meshgrid(xaxis, yaxis, zaxis );
D3D = struct('x',x, 'y',y, 'z',z, 'd', d3d, 'zlab', zlab);
return