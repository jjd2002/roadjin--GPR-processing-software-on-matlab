function OPD = do_modelling()
%
% Callback function to drive the routine "splitstep2dmodel.m"
% 2-D modelling of GPR data with the method of Bitri, A. and Grandjean, G.,
% 1998, "Frequency - wavenumber modelling and migration of 2D GPR data in
% moderately heterogeneous dispersive media", Geophysical Prospecting, 46,
% 287-301.
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%

OPD = discardprocdata;                           % discard the current OPD
% run model ...
[OPD.d,OPD.x,OPD.tt2w,OPD.ntr,OPD.dx,OPD.ns,OPD.dt,cfreq] = splitstep2dmodel; 
% Check if Operation has been canceled or aborted and issue a message
if isempty(OPD.d),
    disp('SPLITSTEP 2D MODEL > Operation aborted - No O/P data returned!');
    return
end

% assign remaining output parameters
% OPD.dt     = OPD.dt*1.0e9; 
OPD.sigpos = 0; 
OPD.zlab   = 'Traveltime (ns)'; 
OPD.xlab   = 'Scan Axis (m)'; 
OPD.TxRx   = 0; 
OPD.Antenna = ['Synthetic Model - ' num2str(cfreq) ' MHz']; 

% Display data
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab);
set(get(gca,'title'),'string',OPD.Antenna,'fonts',14, 'fontw','bold');

% Update processing history
OPD.history = cellstr('Synthetic Model (Split-step method)'); 

return
