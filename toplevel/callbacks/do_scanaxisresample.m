function OPD = do_scanaxisresample(IPD)
%
% Callback function to drive routine "sacnlineresample.m"
% Resample the Scan Axis (change trace spacing).
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%

OPD = discardprocdata; 
% Trap common errors
if isempty(IPD.d), 
   erh = errordlg('No data to process!', 'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end
if isempty(IPD.dx), 
   erh = errordlg('Operation possible only for equally-spaced data!', ...
    'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end
% Proceed 
OPD = IPD; 
[OPD.d,OPD.dx,OPD.ntr,OPD.x] = scanlineresample(IPD.d,IPD.dx,IPD.x);  
if isempty(OPD.d), 
    disp('SCAN AXIS RESAMPLE > Operation aborted - No data returned!');
    return
end
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab); 
% Now resample marker trace id's, albeit roughly
if ~isempty(OPD.markertr),    
    OPD.markertr(:,1) = fix(OPD.markertr(:,1)*OPD.ntr/IPD.ntr) ;
    OPD.markertr(1,1) = IPD.markertr(1,1); 
end;     
% Must resample trace coordinates!
% Case 1. No marker trace info (e.g. Pulse Ekko data, or imported from SU
% or SEGY file formats).
if ~isempty(OPD.xyz.Tx) && isempty(OPD.markertr),
    % Resample Source Antenna coordinates
    dn = IPD.ntr/OPD.ntr;
    x1 = interp1([1:1:IPD.ntr]',IPD.xyz.Tx(:,1),[dn:dn:dn*OPD.ntr]');
    y1 = interp1([1:1:IPD.ntr]',IPD.xyz.Tx(:,2),[dn:dn:dn*OPD.ntr]');
    z1 = interp1([1:1:IPD.ntr]',IPD.xyz.Tx(:,3),[dn:dn:dn*OPD.ntr]');
    OPD.xyz.Tx = [x1(:) y1(:) z1(:)];
    % Resample Receiver Antenna coordinates
    x1 = interp1([1:1:IPD.ntr]',IPD.xyz.Rx(:,1),[dn:dn:dn*OPD.ntr]');
    y1 = interp1([1:1:IPD.ntr]',IPD.xyz.Rx(:,2),[dn:dn:dn*OPD.ntr]');
    z1 = interp1([1:1:IPD.ntr]',IPD.xyz.Rx(:,3),[dn:dn:dn*OPD.ntr]');
    OPD.xyz.Rx = [x1(:) y1(:) z1(:)];
end;     
% Case 2. Monostatic data imported from radar file with marker trace info! 
if ~isempty(OPD.xyz.Tx) && ~isempty(OPD.markertr),
    OPD.xyz.Tx = interpxyz(OPD.dx,OPD.ntr,OPD.markertr); 
    OPD.xyz.Rx = OPD.xyz.Tx; 
end;     
% Processing History
iss = size(OPD.history,1); 
text = ['Resampled Scan Axis, from ' num2str(IPD.ntr) ...
    ' to ' num2str(OPD.ntr) ' traces.'];
OPD.history(iss+1,1) = cellstr(text); 
return
