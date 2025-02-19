function OPD = do_FDTD_modelling()
%
% Callback function to drive the routines "getemproperties.m" and 
% "run_FDTD_model.m" for 2-D, TM-mode, FDTD modeling of reflection
% ground-penetrating radar, implementing the method of James Inrving and
% Rosemary Knight (2006, Computers and Geosciences, 32 1247–1258). 
%
% Copyright (C) 2008, Andreas Tzanis. All rights reserved.
%

OPD = discardprocdata;                           % discard the current OPD

% Import a model and extract its EM property structure
[sig,ep,mu,x,z,t,Fc] = getemproperties;
if isempty(sig),
    disp('FDTD MODELLING > Operation aborted - No O/P data returned!');
    return
end

% Run the simulation
[OPD.d,OPD.x,OPD.tt2w,OPD.dx,OPD.dt] = run_FDTD_model(sig,ep,mu,x,z,t,Fc);
% Check if Operation has been canceled or aborted and issue a message
if isempty(OPD.d),
    disp('FDTD MODELLING > Operation aborted - No O/P data returned!');
    delete(findobj('tag','emsfig'));
    return
end

delete(findobj('tag','emsfig'));
[OPD.ns, OPD.ntr] = size(OPD.d);
% assign remaining aoutput parameters
OPD.dt     = OPD.dt*1.0e9; 
OPD.tt2w   = OPD.tt2w*1.0e9;
OPD.sigpos = 0; 
OPD.zlab   = 'Traveltime (ns)'; 
OPD.xlab   = 'Scan Axis (m)'; 
OPD.TxRx   = 0; 
OPD.Antenna = ['Synthetic Model - ' num2str(Fc*1.0e-6) ' MHz']; 

% Display data
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab);
set(get(gca,'title'),'string',OPD.Antenna,'fonts',14, 'fontw','bold');

% Update processing history
OPD.history = cellstr('Synthetic Model (FDTD simulation)'); 

return
