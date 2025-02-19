function [p, e, attn] = getattenuation(d, tt2w, do_plot, do_model)
%
% GETATTENUATION : * Computes the analytic signal for all traces in the GPR
%                    section, hence their instantaneous power. 
%                  * Computes the median and a mean attenuation function,
%                    i.e. respectively the median and mean instantaneous
%                    power of all traces in the section.
%                  * Computes best fitting models for power-law and
%                    exponential attenuation based on the median
%                    attenuation function.
%                  * Optionally displays the attenuation functions and the
%                    best fitting power-law and exponential decay models.
%
%          Usage : [p, e, attn] = getattenuation(d,tt2w,do_plot,do_model)
%
%         Inputs : 
%              d : [ns x ntr] data array (GPR section)
%           tt2w : [1 x ns] vector of the 2-way traveltime
%        do_plot : Keyword.
%                  = 'no_plot' : Do not display the attenuation curves and
%                                best fitting models
%                  = anything else, or absent : Display the attenuation
%                                curves and best fitting models 
%       do_model : Keyword.
%                  = 'no_tpow' : Do not compute the best fitting power-law
%                                attenuation model
%                  = 'no_expow': Do not compute the best fitting
%                                exponential attenuation model
%                  = anything else, or absent : Compute the best fitting
%                                attenuation models 
%
%        Outputs : 
%              p : [1 x 2] vector, the best fitting power-law attenuation
%                  model. Defaults to [NaN NaN] if do_model = 'no_tpow'.
%              e : [1 x 2] vector, the best fitting exponential attenuation
%                  model. Defaults to [NaN NaN] if do_model = 'no_expow' 
%           attn : [ns x 2] vector holding the median attenuation (1st
%                           column) and the mean attenuation (2nd column).
%
%         Author : Andreas Tzanis,
%                  Department of Geophysics and Geothermy, 
%                  University of Athens
%                  atzanis@geol.uoa.gr
%                 (C) 2005, Andreas Tzanis, all rights reserved
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

if isempty(d), 
    erh = errordlg('No data - the array is empty!','GETATTENUATION: ERROR');
     uiwait(erh);   clear erh; 
     p = [];
     e = [];
     return
end

% Work arrays for fitting attenuation models, global within the scope of
% "getattenuation"
global tt decaycurve

if nargin < 4,
    do_model = 'get_models';
end

if nargin < 3,
    do_plot = 'draw_curves';
end

[ns, ntr] = size(d);
% compute the analytic signal  %%%%%%%%%%%
ns3 = 3*ns;           % this will zero padd traces by 3-fold
ftr = fft(d,ns3);
h   = zeros(ns3,1);
if ns3 >0 && 2*fix(ns3/2)==ns3      % even
    h([1 ns3/2+1]) = 1;
    h(2:ns3/2) = 2;
elseif ns3 > 0                     % odd 
    h(1) = 1;
    h(2:(ns3+1)/2) = 2;
end
ht = ifft( ftr .* (h * ones(1,ntr) ));
ht = abs( ht(1:ns,:) ).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
attn = [median(   ht' )' mean(  ht' )'];
%median_attenuation = smoothh(median_attenuation,10);

decaycurve = attn(:,1)';
tt   = tt2w;
if strcmp(do_model,'no_tpow'),
    p = [NaN NaN];
    tpower = NaN*ones(1,length(tt2w));
else
% Get best fitting power-law attenuation model
    a0   = [max(decaycurve) -2];
    p    = fminsearch(@tpow,a0);
    tpower = p(1)*tt2w.^p(2);
end
if strcmp(do_model,'no_expow'),
    e = [NaN NaN];
    expower = NaN*ones(1,length(tt2w));
else
% Get best fitting exponential attenuation model
    a0   = [max(decaycurve) -1];
    e    = fminsearch(@expow,a0);
    expower = e(1)*exp(e(2)*tt2w);
end

% Check if no plot was requested
if nargin > 2,
    if strcmp(do_plot,'no_plot')
        clear global tt decaycurve
        return
    end
end

% Display the results
attnfig = findobj('tag','attenuationfigure');
if ~isempty(attnfig),
    figure(attnfig);
    clf;
    figuretools;
    matgprwindows('updateadd');
else
   attnfig = figure('name','ATTENUATION CHARACTERISTICS',...
       'tag','attenuationfigure', ... 
       'numbertitle','off', ... 
       'menubar','none', ...
       'CreateFcn',['figuretools; matgprwindows(''setup''); ' ...
                    'matgprwindows(''updateadd'');'],...
       'DeleteFcn',['clear global tt decaycurve; ' ... 
                    'matgprwindows(''updateremove'', ' ... 
                    'findobj(''tag'',''attenuationfigure''));']);
end
semilogy(tt2w, attn(:,1), tt2w, attn(:,2), ...
    tt2w, tpower, tt2w, expower, 'r-.', 'linewidth', 1 )
xlabel('Time (ns)');  ylabel('Instantaneous Power');
grid;
lines = get(gca,'children');
if ~isempty(p),
    pr = num2str(p(2),4);
    plabel = ['Best fitting power-law decay curve ~ t^(' pr ')' ];
else
    plabel = 'Best fitting power-law decay curve: NA' ;
end
if ~isempty(e),
    er = num2str(e(2),3);
    elabel = ['Best fitting exponential decay curve ~ e^(' er 't)'];
else
    elabel = 'Best fitting exponential decay curve: NA';
end
lgd = legend(lines,...
    texlabel(elabel), ...
    texlabel(plabel),...
    'Mean Attenuation', ...
    'Median Attenuation');
set(lgd,'fontsize',10)

attnaxs = uicontrol('Style', 'Popup','tag','attnaxisscale', ... 
    'Parent',findobj('tag','attenuationfigure'), ...
    'String', 'Lin-Lin|Lin-Log|Log-Lin|Log-Log', 'value', 2, ...
    'units','normalized','Position', [0.70 0.04 0.2 0.02], ...
    'tooltip','Select scaling of time and amplitude axes', ...
    'Callback', @do_axsc);

uicontrol('Style', 'pushbutton','tag','attngridonoff', ... 
    'Parent',findobj('tag','attenuationfigure'), ...
    'String', 'Grid', 'SelectionHighlight', 'on', ...
    'units','normalized','Position', [0.132 0.015 0.1 0.05], ...
    'tooltip','Toggle grid on / off', 'Callback', ...
    ['grdstatus = get(gca,''xgrid''); ' ...
     'if strcmp(grdstatus,''on''), ' ...
     '   set(gca,''xgrid'',''off'',''ygrid'',''off''); ' ...
     'end; ' ...
     'if strcmp(grdstatus,''off''), ' ...
     '   set(gca,''xgrid'',''on'',''ygrid'',''on''); ' ...
     'end; ' ]);

return
%
function do_axsc(attnaxs,eventdata)
% Function handle callback to control how to show the spectra of which
% input data trace
axsc  = get(findobj('tag','attnaxisscale'),'value'); 
if axsc == 1, 
     set(gca,'xscale','lin','yscale','lin');
elseif axsc == 2, 
     set(gca,'xscale','lin','yscale','log');
elseif axsc == 3, 
     set(gca,'xscale','log','yscale','lin');
elseif axsc == 4, 
     set(gca,'xscale','log','yscale','log');
end
return

function er = tpow(a)
global decaycurve tt
ll =max(size(decaycurve));
res = a(1)*tt(2:ll).^a(2);
er = norm( log10(decaycurve(2:ll)) - log10(res) );
return

function er = expow(a)
global decaycurve tt
ll =max(size(decaycurve));
res = a(1)*exp(a(2)*tt(2:ll));
er = norm( log10(decaycurve(2:ll)) - log10(res) );
return