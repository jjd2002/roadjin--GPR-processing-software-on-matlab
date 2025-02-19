function saveD3D(D3D) 
%
% SAVED3D : Save the 3-D data volume generated by "Make3DData.m", using
%           MATGPR's generic M3D binary format. The M3D format is  follows:
%--------------------------------------------------------------------------
%     Field             Type       Description 
%  or Variable
%--------------------------------------------------------------------------
%     nx                int16      Size of X-axis (X-dimension)
%     x1                float      Min(X) - Location of first trace
%     dx                float      Trace spacing      
%     ny                int16      Size of Y-axis (Y-dimension)
%     y                 float      Y-axis (vector)    
%     nz                int16      Size of Z-axis (Z-dimension)
%     z1                float      Min(Z) - Minimum value of vertical
%                                  axis. Can be Time or Depth
%     dz                float      Z-spacing (can be Time or Depth)
%     d                 float      [ny x nx x nz] data volume is exported  
%                                  by writting ny slices of size [nx x nz].
%     ls                int16      Size of Z-axis label
%     zlab              uchar      Z-axis label (time or depth).
%
% Usage: saveD3D(D3D)
%
% Input  : D3D, the 3-D data volume generated by "Make3DData.m"
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
% Trap a common error
if isempty(D3D), 
    erh = errordlg('No data to save', 'MATGPR: ERROR'); 
    uiwait(erh);  
    return; 
end; 
% Get file name
[xfname, xpname]= uiputfile( ...
    {'*.m3d','MATGPR M3D-Files (*.m3d)' }, ...
     'Save 3-D Data as ...',ENVAR.currentworkdir); 
if xfname == 0,                               % Cancelled
    return
end
% Export D3D ...
fid = fopen([xpname xfname],'w',ENVAR.endian);
dx = mean(diff(D3D.x(1,:,1)));                % Trace spacing
dz = mean(diff(D3D.z(1,1,:)));                % Depth spacing
nx = size(D3D.x(1,:,1),2);                    % Data dimensions
ny = size(D3D.y(:,1,1),1);
nz = size(D3D.z(1,1,:),3);
fwrite(fid,nx,'int16');                       % Export dim of X-axis     
fwrite(fid,D3D.x(1,1,1),'float');             % Export first trace location
fwrite(fid,dx,'float');                       % Export trace spacing
fwrite(fid,ny,'int16');                       % Export dim of Y-axis     
fwrite(fid,D3D.y(:,1,1),'float');             % Export Y-axis
fwrite(fid,nz,'int16');                       % Export dim of Z-axis     
fwrite(fid,D3D.z(1,1,1),'float');             % Min Z-value                
fwrite(fid,dz,'float');                       % Export Z-spacing
for i=1:ny                                    % Export 3-D data volume
    fwrite(fid,D3D.d(i,:,:),'float');
end
ls = max(size(D3D.zlab));                     % Export dim of Z-axis label
fwrite(fid,ls,'int16');
fwrite(fid,D3D.zlab,'uchar');                 % Export Z-axis label
fclose(fid);
return
