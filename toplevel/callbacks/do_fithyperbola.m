function do_fithyperbola(IPD)
%
% Callback function to drive the routine "fitdiffractionhyperbola.m"
% Interactively model observed diffraction front hyperbolae 
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%

global VS vhalfspace 
% Trap common errors ...
if isempty(IPD.d),                               
   erh = errordlg('No data to process!', 'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end
if isempty(IPD.dx),                              
    erh = errordlg('Operation possible only for equally-spaced data!', ...
        'MATGPR : ERROR');
    uiwait(erh); 
    return;
end;
% Run program 
fitdiffractionhyperbola(IPD); 
waitfor(findobj('Tag','vel_handle'));
if ~isempty(vhalfspace),
    VS.v1d = [vhalfspace 0]; 
else
    VS.v1d = [];
end;
showinfo(IPD); 
clear global vhalfspace 
return
