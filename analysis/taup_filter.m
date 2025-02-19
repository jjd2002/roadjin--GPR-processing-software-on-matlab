function [drf, tp, p] = taup_filter(d, t, dt, x, dx)
%
% TAUP_FILTER : Model and De-noise data by means of the linear Tau-P
%               transform. The options are zone-pass or zone-stop (mute)
%               filtering. This works like a band pass or stop filter. A
%               polygonal zone is defined a-priori, or graphically by
%               pointing and clicking at the desired positions on the Tau-P
%               panel. "Zone pass" means that ONLY the Tau-P contributions
%               IN and ON the sides of the zone will be retained. "Zone
%               stop" means that ONLY the Tau-P contributions OUTSIDE the
%               zone will be retained. 
%            => A necessary condition is that traces are equally spaced.
%            => THIS PROGRAM IS PART OF THE MATGPR SUITE.

%     Usage : [drf, tp, p] = taup_filter(d, t, dt, h, dx) 
%   
%    Inputs : 
%         d : The 2-D input data matrix (radargram)
%         t : Two-way traveltime vector in ns
%        dt : Sampling rate (in ns)
%         x : The scan-axis vector (trace locations along survey line) in m
%        dx : Sampling rate along rows (trace spacing in meters)
%
%    Output : 
%       drf : Data model after inverse Tau-P transform
%        tp : The Tau-P panel
%         p : The vector of slopes P (slowness)
%
%  Requires : fwd_taup_tx, inv_taup_tx, setpolygon, isinpoly
%
% Author    : Andreas Tzanis,
%             Department of Geophysics, 
%             University of Athens
%             atzanis@geol.uoa.gr
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

% Get environmental variables
global ENVAR 
global rdtype 
%%%%%    Inquire filter type
rdtype = menu('SELECT FILTER TYPE','Zone Pass',...
    'Zone Stop', 'Cancel');
if rdtype == 3,
    drf = [];  tp=[]; p=[]; 
    return
end
% Determine whether the pass / stop zone coordinates are predefined (import
% from file), or will be designed on-screen
inmode = menu('PLEASE DECIDE INPUT METHOD',...
    'Import Zone Coordinates from File','Design Filter on Screen',...
    'Cancel');
if inmode == 3,
    drf = []; tp=[]; p=[];
    return
end

%%%%%   Useful parameters
d_norm = norm(d);                                  % to be used for scaling
fn = 1./(2*dt);                                  % Nyquist frequency in GHz
fmax = 0.75*fn;          % max frequency of interest as fraction of Nyquist
pmax = 1/(2*dx*fmax);                                           % max slope
pmin = -pmax;
dp = 1/(max(x)*fmax);
p    = pmin:dp:pmax;

% Forward transform t-x --> tau-p
tp = fwd_taup_tx(d, dt, x, dx, p);

% Design pass/ stop zone
%%%%%%% Use adaptive figure sizing and posistioning %%%%%%%%%%%%%%%%%%%%%%%
scrsz  = get(0,'screensize');
ffkpos = round([440*scrsz(3)/1680 200*scrsz(4)/1050 ...
    800*scrsz(3)/1680 700*scrsz(4)/1050]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ffk = figure('Numbertitle','off','tag','taupfiltfig',...
    'Name','Tau-P Panel','menubar','none', ...
    'position', ffkpos);
pcolor(p, t, tp);
shading('flat');
set(gca,'ydir','reverse');
if ~isempty(ENVAR) && ~isempty(ENVAR.colormap),
    colormap(ENVAR.colormap);
end
imagecolors;
ylabel('Tau (ns)');
xlabel('slowness (p, ns/m)')
set(gca,'nextplot','add')
%%% Zone data imported from file
if inmode == 1,
    [rdname, rdpath]= uigetfile('*.dat;*.txt',...
        'Please give file with coordinates of Pass/Stop Zone', ...
        ENVAR.currentworkdir);
    if rdname == 0,                           % Canceled
        delete(ffk);
        drf = []; tp=[]; p=[]; 
        return;
    end;
    poly = importdata([rdpath rdname]);   % get p or q and tau vertices of zonal polygon
    slope    = poly(:,1);
    tau      = poly(:,2);
    clear poly
    %
    patch3 = patch(slope,tau,'k');
    set(patch3,'facealpha',0.1)
    reply = questdlg('Is this OK ? ');
    if strcmp(reply,'Cancel') || strcmp(reply,'No'),
        delete(ffk);
        drf = []; tp=[]; p=[]; 
        return
    end

    %%% Zone will be designed on screen
elseif inmode == 2,
    % Tell user how to enter the pass/stop zones using mouse
    txt = cell(1);
    txt(1) = cellstr('To set vertices of polygonal pass / stop zone,');
    txt(2) = cellstr('click LEFT mouse button. To correct mistakes ');
    txt(3) = cellstr('click the MIDDLE button. When finished, click ');
    txt(4) = cellstr('the RIGHT button.');
    msg = msgbox(txt);
    pos=get(msg,'position');
    set(msg,'pos',[20 20 pos(3) pos(4)]);
    pause(0.5)
    reply = 'No';
    while strcmp(reply,'Yes')~=1,
        figure(ffk)
        % k, f are the ray parameter and time coordinates of the polygon
        % vertices respectively
        [slope, tau] = setpolygon(' ns/m', ' ns');
        if isempty(slope) || isempty(tau),
            delete(ffk);
            if exist('msg','var') && ishandle(msg),
                delete(msg)
            end
            drf = []; tp=[]; p=[];
            return
        end
        patch3 = patch(slope,tau,'k');
        set(patch3,'facealpha',0.1)
        %%%%%   Make sure that pass-zone has been set to satisfaction
        reply = questdlg('Is this OK ? ');
        if strcmp(reply,'Cancel')==1,
            delete(ffk);
            if exist('msg','var') && ishandle(msg),
                delete(msg)
            end
            drf = []; tp=[]; p=[]; 
            return
        end
        if strcmp(reply,'No')==1,
            delete(patch3);
            
            pause(0.5)
        end
    end
end
%%% Clear first message to make room for the second!
if exist('msg','var') && ishandle(msg),
    delete(msg)
end
 
%%% Set up ray parameter and time grids for "inpolygon" and "isinpoly"
[P, TAU] = meshgrid(p, t);
%%% The MATLAB function "inpolygon" may not exist for MATLAB releases
%%% earlier that 13! Check and if not, switch to Kirill Pankratov's
%%% "isinpoly".  
find_inpolygon = which('inpolygon');
if ~isempty(find_inpolygon),
    IN  = inpolygon(P,TAU,slope,tau);               
elseif isempty(find_inpolygon),
    disp('TAUP_TX > INPOLYGON not found! Working with ISINPOLY');
    IN  = inpoly(P,TAU,slope,tau);
end
OUT = ~IN;                                     
if rdtype==1,                                                   % zone pass 
    tpf = tp.*IN + 0.01*tp.*OUT;
elseif rdtype==2,
    tpf = 0.01*tp.*IN + tp.*OUT;                                % zone stop 
end
%clear P TAU IN OUT             
delete(ffk); clear ffk          
    
% Inverse transform tau-p/ tau-q --> t-x
drf = inv_taup_tx(tpf, dt, x, dx, p);

% Scaling
m_norm = norm(drf);                                            % model norm
drf = d_norm*(drf/m_norm);                            % scale model to data

% Display results and decide what to do...
ffk = figure('Numbertitle','off','tag','rdfiltfig',...
    'Name','MODEL AND RESIDUALS','menubar','none', ...
    'position', ffkpos);
figuretools;
imagecolors;
% Create "actions" menu
uimenu('Label',' ACTIONS', 'Tag', 'ffkmenu');
uimenu(findobj('Tag','ffkmenu'), 'tag', 'tpunscalemodel', ...
    'Label', 'Unscale Model', ... 
    'Callback', {@unscalemodel, d_norm, m_norm} );
uimenu(findobj('Tag','ffkmenu'), 'tag', 'tpscalemodel', ...
    'Label', 'Re-scale Model', ... 
    'Enable', 'off', ...
    'Callback', {@rescalemodel, d_norm, m_norm} );
uimenu(findobj('Tag','ffkmenu'), ...
    'Label', 'Keep Model', ... 
    'Callback', @keepmodel );
uimenu(findobj('Tag','ffkmenu'), ...
    'Label', 'Keep Residuals', ... 
    'Callback', @keepresiduals ); 
uimenu(findobj('Tag','ffkmenu'), ...
    'Label', 'Discard and return', ... 
    'Callback', @discardall ); 
% Display data model
cax = get(findobj('tag','currentdataaxes'),'clim');
subplot(211);
imagesc(x, t, drf);
set(gca,'clim',cax);
ylabel('t (ns)');  xlabel('x (m)');  title('Model');
% Display residuals
subplot(212);
imagesc(x, t, d-drf);
set(gca,'clim',cax);
ylabel('t (ns)');  xlabel('x (m)');  title('Residuals');
if ~isempty(ENVAR) && ~isempty(ENVAR.colormap),
    colormap(ENVAR.colormap);
end
waitfor(ffk)
% Return result
%drf = result;

return

function unscalemodel(obj, evt, dnorm, mnorm)
x   = evalin('caller','x');             % get variables from caller routine
t   = evalin('caller','t');
d   = evalin('caller','d');
drf = evalin('caller','drf');
drf = mnorm*(drf/dnorm);                                    % unscale model
assignin('caller','drf', drf);            % replace model in caller routine
figure(findobj('tag','rdfiltfig'));
cax = get(findobj('tag','currentdataaxes'),'clim');
% Display model
subplot(211);
imagesc(x, t, drf);    set(gca,'clim',cax);
ylabel('t (ns)');  xlabel('x (m)');  title('Model');
% Display residuals
subplot(212);
imagesc(x, t, d-drf);    set(gca,'clim',cax);
ylabel('t (ns)');  xlabel('x (m)');  title('Residuals');
set(findobj('tag', 'tpunscalemodel'),'enable','off')
set(findobj('tag', 'tpscalemodel'),'enable','on')
return

function rescalemodel(obj, evt, dnorm, mnorm)
x   = evalin('caller','x');             % get variables from caller routine
t   = evalin('caller','t');
d   = evalin('caller','d');
drf = evalin('caller','drf');
drf = dnorm*(drf/mnorm);                                    % unscale model
assignin('caller','drf', drf);            % replace model in caller routine
figure(findobj('tag','rdfiltfig'));
cax = get(findobj('tag','currentdataaxes'),'clim');
% Display model
subplot(211);
imagesc(x, t, drf);    set(gca,'clim',cax);
ylabel('t (ns)');  xlabel('x (m)');  title('Model');
% Display residuals
subplot(212);
imagesc(x, t, d-drf);    set(gca,'clim',cax);
ylabel('t (ns)');  xlabel('x (m)');  title('Residuals');
set(findobj('tag', 'tpunscalemodel'),'enable','on')
set(findobj('tag', 'tpscalemodel'),'enable','off')
return

function keepmodel(obj, evt)
global rdtype
if rdtype == 1,
    rdtype = 11;
end
if rdtype == 2.
    rdtype = 21;
end
delete(findobj('tag','rdfiltfig'))
return

function keepresiduals(obj, evt)
global rdtype
d = evalin('caller','d');
m = evalin('caller','drf');
m = d - m; 
assignin('caller','drf', m);
if rdtype == 1,
    rdtype = 12;
end
if rdtype == 2.
    rdtype = 22;
end
delete(findobj('tag','rdfiltfig'))
return

function discardall(obj, evt)
global rdtype
m = evalin('caller','drf');
m = []; 
assignin('caller','drf', m);
rdtype = 3;
delete(findobj('tag','rdfiltfig'))
return
