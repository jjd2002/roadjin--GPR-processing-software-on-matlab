function OPD = do_fx_decon(IPD)
%
% Callback function to drive the routine "fx_decon.m"
% Perform F-X Deconvolution 
%
% Copyright (C) 2008, Andreas Tzanis. All rights reserved.
%

OPD = discardprocdata; 
% Trap common errors
if isempty(IPD.d), 
   erh = errordlg('No data to process!', 'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end

% Ask how deconvolution parameters will be given. If the data is not
% equally spaced, do not bother ...
if ~isempty(IPD.dx),
    ask = questdlg('Filter parameters given in # samples or in metres?',...
        'Prediction Operator parameters', ...
        'Samples','metres','Cancel','Samples');
    if strcmp(ask,'Cancel'),
        disp('F-X DECONVOLUTION > Operation aborted - No data returned!');
        return
    end
end

% Get filter parameters
% lf,     Length of the filter
% mu,     Prewhitening in %  
% flow,   min  freq. in the data in Hz
% fhigh,  max  freq. in the data in Hz
fhigh = 1./(2*IPD.dt);             % Nyquist frequency
nf    = 2^nextpow2(IPD.ns);
flow  = fhigh/(nf-1);
defaultans = {'', '', num2str(flow), num2str(fhigh)};
answer = inputdlg({'Prediction Operator Length' ... 
    '% Prewhitening' ...
    'Minimum frequency in data (GHz)' ...
    'Maximum frequency in data (GHz)' },...
    'Decovolution parameters',1,defaultans);
if isempty(answer),
    disp('F-X DECONVOLUTION > Operation aborted - No data returned!');
    return
end
ansr  = checkcomma(answer); 
lf    = str2num(ansr(1,:));
mu    = str2num(ansr(2,:));
flow  = str2num(ansr(3,:));
fhigh = str2num(ansr(4,:));
if ~isempty(IPD.dx) && strcmp(ask,'metres'),
    lf = floor(lf/IPD.dx);
end
% More error trapping
if lf >= IPD.ntr/2, 
    erh = errordlg(['Operator cannot be longer than ' ...
                    'half the scan axis length!          '], ...
                'F-X DECONVOLUTION : ERROR');
    uiwait(erh);
    return
end

% Now proceed ...
OPD = IPD; 
OPD.d = fx_decon(IPD.d,IPD.dt,lf,mu,flow,fhigh);

% Display result
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab); 
% Update history
iss = size(OPD.history,1); 
if ~isempty(IPD.dx),
    text = ['F-X Deconvolution: Operator length = ' num2str(lf*IPD.dx,4) ' m.'];
else
    text = ['F-X Deconvolution: Operator length = ' num2str(lf) ' traces.'];
end
OPD.history(iss+1,1) = cellstr(text); 
title(text)
return

