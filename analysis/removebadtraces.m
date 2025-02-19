function [dout, bad] = removebadtraces(x,t,d)
%
% REMOVEBADTRACES : Pick bad traces by pointing and clicking and replace
% them with the interpolants computed from their near neighbours. Note,
% however, that:
% 1. Because it uses interpolation, this function works best for isolated
%    bad traces, or small clusters adjacent bad traces. In the latter case
%    the clusters should better be separated by extended groups of
%    consecutive good traces.
% 2. For the same reason, it is NOT RECOMMENDED for removing extended
%    groups of bad traces, as it will probably fail to produce reliable
%    interpolants. 
%
% Usage : [dout, bad] = removebadtraces(x,t,d);
%         dout = removebadtraces(x,t,d);
%
% Input :
%     d : The [ns x ntr] GPR data section        
%     t : [1 x ns] vector of time coordinates (2 way traveltime)
%     x : [1 x ntr] vector of scan axis coordinates
%
%Output :
%  dout : The [ns x ntr] GAPR section without bad traces
%   bad : Vector holding the sequential numbers of bad traces
%
%Author : Andreas Tzanis
%         Department of Geophysics, 
%         University of Athens
%         atzanis@geol.uoa.gr
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
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
[ns, ntr] = size(d);

%%% Create badtraces figure;
badtrfig = figure('name','REMOVE BAD TRACES','tag','badtrfigure',...
    'menubar','none', 'numbertitle','off', ...
    'position',[200   300   800   450], ...
    'CreateFcn',['figuretools; imagecolors; ' ... 
                 'matgprwindows(''setup''); ' ...
                 'matgprwindows(''updateadd'');'],...
    'DeleteFcn',['matgprwindows(''updateremove'', ' ...
                 'findobj(''tag'',''badtrfigure''));']);

% Draw data
if strcmp(lower(ENVAR.DisplayOption),'image')==1,
    imagesc(1:1:ntr,t,d);
    colormap(ENVAR.colormap);

elseif strcmp(lower(ENVAR.DisplayOption),'wiggle') || ...
        strcmp(lower(ENVAR.DisplayOption),'vararea') || ...
        strcmp(lower(ENVAR.DisplayOption),'wig&img') || ...
        strcmp(lower(ENVAR.DisplayOption),'var&img'),

    wiggledisplay(d, 1:1:ntr, t, ENVAR.DisplayOption, ...
        ENVAR.maxwigs, ENVAR.wigglescale);

else
    erh = errordlg('Cannot figure out the display option',...
        'REMOVE BAD TRACES: ERROR');
    uiwait(erh)
    return
end
%%% Set hold on
set(gca,'nextplot','add');
%%% Add axes labels
xlabel('# Traces');  ylabel('Traveltime (ns)');

%%%% Display help message %%%
msgtxt = cell(1);
msgtxt(1) = cellstr('To pick a bad trace click the LEFT button. ');
msgtxt(2) = cellstr('To correct mistakes click the MIDDLE button. ');
msgtxt(3) = cellstr('To finish, click the RIGHT mouse button! ');
msg = msgbox(msgtxt);
pos=get(msg,'position');
set(msg,'pos',[20 20 pos(3) pos(4)]);
pause(0.5)
%%% Now make sure that the right figure is on focus
figure(findobj('tag','badtrfigure'));  

%%% Proceed - pick and mark bad traces      
bad   = [];                       % bad trace numbers will be stored here
but   = 1;                        % mouse button pressed key
nmark = 0;
lin2  = [];
while but==1 || but == 2
    [xp,yp,but] = getpoint(' ', ' ns');
    xp = round(xp);
    if but ==1,
        nmark = nmark+1;
        bad = [bad xp];
        lin2(nmark) = plot([xp xp],[t(1), t(ns)], 'r-','linewidth',2);
    end
    if but == 2,                 % backstep to correct mistakes
        if nmark <= 0,
            continue
        elseif nmark == 1,
            bad = [];
            delete(lin2(nmark))
            nmark = nmark -1;
        else
            bad = bad(1:nmark-1);
            delete(lin2(nmark))
            nmark = nmark -1;
        end
    end
end                              % while loop

%%% Kill help message
if exist('msg') && ishandle(msg),              
    delete(msg);
end

%%% Process (replace) bad traces
bad = sort(bad);
notbad = setxor(1:1:ntr,bad);
%%% Interpolate 
good = zeros(ns, ntr);
for i=1:ns
    good(i,:) = interp1(x(notbad), d(i,notbad), x, 'cubic');
end

%%% Replace bad traces
dout        = d;
dout(:,bad) = good(:,bad);

%%% Remove figure and exit
delete(badtrfig);

return
