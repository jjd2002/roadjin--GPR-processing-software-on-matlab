function [dout, P] = KarhunenLoeve(d, x, t)
%
% KARHUNENLOEVE : Karhunen Loeve filtering to enhance the lateral
%                 coherence of GPR events. Uses "svds" to compute the
%                 P largest singular values and associated singular
%                 vectors of the input data. The KL transfrom is re-
%                 synthesized from the largest P singular values /vectors
%                 and is a low-dimensional (smooth) approximation of the
%                 input data. 
%                 ==> The routine displays the reconstructed data and the
%                 residuals after subtracting the reconstructed data from
%                 the input data. By means of a specialized (ACTIONS) menu,
%                 the user may decide whether to keep the model, or the
%                 residuals as the output of the K-L filtering operation.
%
%         Usage : [dout, P] = KarhunenLoeve(d, x, t)
%
%         Input : 
%             d : The [ns X ntr] input (full dimensional) data matrix 
%             x : The [1 x ntr ] scan axis vector (trace locations).
%             t : The [1 x ns) 2-way traveltime vector
%
%        Output : 
%          dout : => The low-dimensional approximation to the input data "d"
%                 (smoothed GPR section), or,
%                 => The residuals after subtracting the reconstructed data
%                 from the input data
%             P : Size of the approximating low-dimensional sub-space.
%
%      Requires : checkcomma.m
%
%        Author : Andreas Tzanis,
%                 Department of Geophysics, 
%                 University of Athens
%                 atzanis@geol.uoa.gr
%
% Copyright (C) 2005,2008 Andreas Tzanis. All rights reserved.
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

% Get # of eigenvector to be used
answer = inputdlg(['变换中使用的最大特征向量的数量'], 'Karhunen - Loeve filtering',1);
if isempty(answer),
    dout = [];  P=[];
    return
end
ansr = checkcomma(answer); 
P = str2num(ansr(1,:));
% Compute the first P singular values / vectors using "svds" 
disp('KARHUNEN LOEVE > Please wait while computing the SVD')
msh = msgbox('Please wait while computing the SVD and synthesizing output',...
    'KarhunenLoeve','help');
[U,S,V] = svds(d,P);
% Synthesize dout from the first P singular values / vectors
dout    = U*S*V';
%%% let message window go
if ishandle(msh),
    close(msh);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Plot model and residuals and decide how to proceed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Adaptive figure sizing and posistioning 
scrsz  = get(0,'screensize');
klfpos = round([440*scrsz(3)/1680 200*scrsz(4)/1050 ...
    800*scrsz(3)/1680 700*scrsz(4)/1050]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
klf = figure('Numbertitle','off','tag','klffiltfig',...
    'Name','KL FILTER: RECONSTRUCTION AND RESIDUALS','menubar','none', ...
    'position', klfpos);
figuretools;
imagecolors;
% Create "actions" menu
uimenu('Label',' ACTIONS', 'Tag', 'klfmenu');
uimenu(findobj('Tag','klfmenu'), ...
    'Label', 'Keep Reconstruction', ... 
    'Callback', 'delete(findobj(''tag'',''klffiltfig''));' );
uimenu(findobj('Tag','klfmenu'), ...
    'Label', 'Keep Residuals', ... 
    'Callback', @keepresiduals ); 
uimenu(findobj('Tag','klfmenu'), ...
    'Label', 'Discard and return', ... 
    'Callback', @discardall ); 
% Display data model
cax = get(findobj('tag','currentdataaxes'),'clim');
subplot(211);
imagesc(x, t, dout);
set(gca,'clim',cax);
ylabel('t (ns)');  xlabel('x (m)');  title('RECONSTRUCTED DATA');
% Display residuals
subplot(212);
imagesc(x, t, d-dout);
set(gca,'clim',cax);
ylabel('t (ns)');  xlabel('x (m)');  title('RESIDUALS');
if ~isempty(ENVAR) && ~isempty(ENVAR.colormap),
    colormap(ENVAR.colormap);
end
waitfor(klf)
%%% that's all

return

function keepresiduals(obj, evt)
global rdtype
d = evalin('caller','d');
m = evalin('caller','dout');
m = d - m; 
assignin('caller','dout', m);
delete(findobj('tag','klffiltfig'))
return

function discardall(obj, evt)
global rdtype
m = []; 
assignin('caller','dout', m);
delete(findobj('tag','klffiltfig'))
return
