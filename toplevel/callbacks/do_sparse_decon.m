function OPD = do_sparse_decon(IPD)
%
% Callback function to drive the routine "sparse_decon.m"
% Sparse Spike Deconvolution 
%
% Copyright (C) 2009, Andreas Tzanis. All rights reserved.
%

global ENVAR
OPD = discardprocdata; 
% Trap common errors
if isempty(IPD.d), 
   erh = errordlg('No data to process!', 'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end

% Parameters
% mu,        regulatization factor  
% itermax,   maximu # of iterations      
% w,         source wavelet        
defaultans = {'0.1', '20', '1', ''};
answer = inputdlg({'Regularization factor (mu)' ... 
    'Maximum Iterations' ...
    'Wavelet: 1 = Blackman-Harris,  2 = Ricker' ...
    'Antenna Frequency (MHz)' },...
    'Decovolution parameters',1,defaultans);
if isempty(answer),
    disp('SPARSE DECONVOLUTION > Operation aborted - No data returned!');
    return
end
if isempty(answer{4}),
   erh = errordlg('Please provide Antenna Frequency', 'SPARSE DEC0N: ERROR'); 
   uiwait(erh);  
   return
end    
ansr    = checkcomma(answer); 
mu      = str2num(ansr(1,:));
itermax = str2num(ansr(2,:));
wtype   = str2num(ansr(3,:));
freqc   = str2num(ansr(4,:));
if wtype == 1,
    w = blackharrispulse(freqc/1000,IPD.tt2w); 
    ii = w ~= 0; 
    w=w(ii); 
    w=w(:);
    wavelet = 'Blackman-Harris';
elseif wtype == 2,
    w = ricker(freqc/1000,IPD.dt);
    w = w(:);
    wavelet = 'Ricker';
end

% Now proceed ...
result = [];
OPD = IPD; 
[refl, dp] = sparse_decon(IPD.d, w, mu, itermax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Plot reflectivit and model and decide how to proceed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Adaptive figure sizing and posistioning 
scrsz  = get(0,'screensize');
spdcpos = round([440*scrsz(3)/1680 200*scrsz(4)/1050 ...
    800*scrsz(3)/1680 700*scrsz(4)/1050]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spdc = figure('Numbertitle','off','tag','spdcfiltfig',...
    'Name','SPARSE DECONVOLUTION - DATA MODEL AND REFLECTIVITY','menubar','none', ...
    'position', spdcpos);
figuretools;
imagecolors;
% Create "actions" menu
uimenu('Label',' ACTIONS', 'Tag', 'spdcmenu');
uimenu(findobj('Tag','spdcmenu'), ...
    'Label', 'Keep Model', ... 
    'Callback', @keepDp );
uimenu(findobj('Tag','spdcmenu'), ...
    'Label', 'Keep Reflectivity', ... 
    'Callback', @keepRefl ); 
uimenu(findobj('Tag','spdcmenu'), ...
    'Label', 'Discard and Return', ... 
    'Callback', @discardall ); 
% Display data model
subplot(211);
imagesc(OPD.x, OPD.tt2w, dp);
ylabel('t (ns)');  xlabel('x (m)');  title('MODEL');
% Display residuals
subplot(212);
imagesc(OPD.x, OPD.tt2w, refl);
ylabel('t (ns)');  xlabel('x (m)');  title('REFLECTIVITY');
if ~isempty(ENVAR) && ~isempty(ENVAR.colormap),
    colormap(ENVAR.colormap);
end
waitfor(spdc)

% Assign result to OPD structure
if strcmp(result,'Reflectivity'),
    OPD.d = refl;
elseif strcmp(result,'Model'),
    OPD.d = dp;
else
    OPD.d = [];
end

% Display result
if isempty(OPD.d), 
    disp('SPARSE DECON > Operation aborted - No O/P data returned!');
    return
end
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab); 
% Update history
iss = size(OPD.history,1); 
text = ['Sparse Deconvolution: ' result ', Wavelet=' wavelet ', mu=' num2str(mu) ];
OPD.history(iss+1,1) = cellstr(text); 
title(text)
return

function keepDp(obj, evt)
assignin('caller','result', 'Model');
delete(findobj('tag','spdcfiltfig'))
return

function keepRefl(obj, evt)
assignin('caller','result', 'Reflectivity');
delete(findobj('tag','spdcfiltfig'))
return

function discardall(obj, evt)
assignin('caller','result', 'Cancel');
delete(findobj('tag','spdcfiltfig'))
return
