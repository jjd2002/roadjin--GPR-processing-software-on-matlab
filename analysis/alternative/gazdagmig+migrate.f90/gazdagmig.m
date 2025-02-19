function dmig  = gazdagmig(d, dt, dx, vofh)
%
%  GAZDAGMIG  : Slow but accurate Phase-shifting migration for constant 
%               or layered velocity halfspaces
%
%  Inputs     : d is the common-offset section
%             : dt is the sampling interval in time
%             : dx is trace spacing. 
%             : vofh(nlay,2) is the 1-D velocity model of "nlay" velocity - 
%               thickness pairs
%  Outputs    : dmig is the migrated section
%
%  Author     : Andreas Tzanis,
%               Dept. of Geophysics,   University of Athens
%               atzanis@geol.uoa.gr
%  Created    : September 2003  
%

[ns,ntr]=size(d);
%%%%%   Check Velocity structure
[nlay, vhpairs] = size(vofh);
if vhpairs ~=2,
    erh = errordlg('Error in the structure of the velocity model',...
        'GAZDAGMIG : ERROR');
    uiwait(erh)
    dmig = [];
    return
end
layer_velocity  = vofh(:,1)';
layer_thickness = vofh(:,2)';
for i=1:nlay
    if layer_velocity(i) <= 0 || layer_velocity(i) > 0.2998,
        errordlg('Impossible velocity value found! Please try again!',...
            'GAZDAGMIG : ERROR')
        dmig = [];   
        return
    end
end
%%%  Begin
tic                                      % start counting processing time
if nlay == 1,                            %%%%%   Case of uniform halfspace
    vmig = layer_velocity(1);
elseif nlay > 1,                         %%%%%   Case of layered halfspace
    tt2w      = 0:dt:(ns-1)*dt;          % two way traveltime
    firstz = 0.0;
    firstt = 0.0;
%%% Compute migration velocity vmig(t) from v(z) 
    zmax = max(vofh(:,1))*tt2w(ns)/2;    % max possible penetration
    delz = (zmax-firstz)/(ns-1);
    z    = firstz : delz : zmax;
    nz        = length(z);
%%% Compute t(z) from V(z)
    temptz(1) = 2.0*firstz/layer_velocity(1);
    for iz=1:nlay-1; 
        temptz(iz+1) = temptz(iz) + 2.0*layer_thickness(iz)/layer_velocity(iz);
    end
    temptz(nlay + 1) = temptz(nlay) + ...
        2.0*(zmax - sum(layer_thickness(1:nlay)))/layer_velocity(nlay);
    tofz  = interp1([0 cumsum(layer_thickness(1:nlay-1)) zmax] ,temptz,z);
%%% Compute z(t) from t(z)
    vfz   = layer_velocity(1);             % initial velocity
    vlz   = layer_velocity(nlay-1);        % final velocity at depth z
    lt    = firstt+(ns-1)*dt;
    lz    = firstz+(nz-1)*delz;
    zoft  = yxtoxy(nz,delz,firstz,tofz,ns,dt,firstt,0.0,0.0);
    ii = find(tt2w < tofz(1));             % take care of out of range values
    if ~isempty(ii),
        zoft(ii) = 0.5*tt2w(ii)*vfz;
    end
    ii=find(tt2w >= tofz(nz));
    if ~isempty(ii),
        zoft(ii) = lz + 0.5*(tt2w(ii) - tofz(nz))*vlz;
    end
%%% Compute vmig(t) from z(t) 
    vmig = [];
    for it = 1:ns-1,
        vmig(it) = 2.0*(zoft(it+1) - zoft(it))/dt;
    end
    vmig(ns) = vmig(ns-1);
end                                         %%%%%   if nlay loop

%%% Fourier transform from time-space (TX) to frequency-wavenumber (FKx)
%%% domain
fk  = fft2(d);      
fk  = fftshift(fk);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute image in the time-wavenumber (TKx) domain using the         
%%% fast fortran program migrate.f90                                   
[imgr,imgi] = external_f90(ns, ntr, dt, dx, nlay, vmig, fk);
if isempty(imgr)
    dmig = [];
    return
end
img = imgr + sqrt(-1)*imgi;                      % complex image         
clear imgr imgi;                                                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Transform from time-wavenumber (TKx) to time-space (TX) domain and 
%%% get migrated image
h = waitbar(0,'Transforming TK_x -> TX');
dmig = zeros(ns,ntr);
for i=1:ns,
    dmig(i,:) = ifft(ifftshift(img(i,:)));
    waitbar(i/ns,h);
end
dmig = real(dmig);
close(h);
time  = toc;
disp(['GAZDAGMIG finished in ' num2str(time) ' seconds'])
return
%
function [imgr,imgi] = external_f90(ns, ntr, dt, dx, nlay, vmig, fk)
%
% Driver for the external Fortran 90 program MIGRATE.EXE, which will  
% do the heavy numerical work (image computations).
%
%  Author     : Andreas Tzanis,
%               Dept. of Geophysics,   University of Athens
 
%%%%% Export FKx data for migrate.f90 
OutDir = tempdir;
fid    = fopen([OutDir 'for_migration.dat'],'w');
if nlay == 1,
    fwrite(fid,'cvmigg','char');             
elseif nlay >1,
    fwrite(fid,'lvmigg','char');             
end
fwrite(fid,ns,'int');      
fwrite(fid,ntr,'int');
fwrite(fid,dt,'float');
fwrite(fid,dx,'float');
if length(vmig) == 1,                  % uniform halfspace
    fwrite(fid,vmig,'float');
elseif length(vmig) == ns,             % layered halfspace
    fwrite(fid,vmig(1:ns),'float');
end
for i=1:ns; 
    fwrite(fid,real(fk(i,:)),'float'); 
end
for i=1:ns; 
    fwrite(fid,imag(fk(i,:)),'float'); 
end
fclose(fid);
%%%%%   Notify user 
txt = cell(1);
txt(1) = cellstr('Working, please wait! The heavy numerical work is being done with');
txt(2) = cellstr('the external program MIGRATE.EXE. Progress can be monitored in the');
txt(3) = cellstr('command window');
msg = msgbox(txt,'GAZDAGMIG: NOTIFICATION','help');
pos = get(msg,'position');
set(msg,'pos',[20 20 pos(3) pos(4)]);

%%%%%   Construct phase-migtated image (in TKx domain)
migratepath = which('migrate.exe');
[busted, result] =...
    dos([migratepath ...
    ' ' OutDir 'for_migration.dat  ' OutDir 'from_migration.dat'],'-echo');
if busted,
    disp('GAZDAGMIG > Ungraceful termination of MIGRATE.EXE')
    disp('GAZDAGMIG > due to system error.')
    imgr = []; imgi = [];
    if exist('msg') && ishandle(msg), %#ok<EXIST>
        delete(msg)
    end
    return
end
%%%%%   Import the image in real and imaginary parts
fopen([OutDir 'from_migration.dat'],'r');
imgr = fread(fid,[ns,ntr],'float');
imgi = fread(fid,[ns,ntr],'float');
fclose(fid);
delete([OutDir 'from_migration.dat']);

%%%%%   Remove the message box
if exist('msg') && ishandle(msg), %#ok<EXIST>
    delete(msg)
end
return
