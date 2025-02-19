function dout = gaininvdecay(d, tt2w)
%
% GAININVDECAY: Applies a gain function which exactly compensates for the
%               mean or median amplitude attenuation of a 2-D GPR section 
%               * Computes the analytic signal for all traces in the GPR
%                 section, hence their instantaneous amplitude. 
%               * Computes a median and a mean amplitude attenuation
%                 function, i.e. respectively the median and mean
%                 instantaneous amplitude of all traces in the section.
%               * Computes best fitting attenuation model with a function
%	                  A =  c(1)*exp(-a(1)*t) + ... + c(n)*exp(-a(n)*t)
%	              with n linear parameters and n nonlinear parameters.
%               * Displays the result and aAllows the user to experiment
%                 with the order n of the fitting function and the median
%                 or mean decay curve.
%               * At the user's request applies gain by multiplying each
%                 trace with the normalized inverse of the amplitude decay
%                 model 1/y/max(y) and exits.
%
%       Usage : dout = gaininvdecay(d, tt2w)
%
%         Inputs : 
%              d : the GPR section
%           tt2w : the 2-way traveltime
%
%        Outputs : 
%           dout : The amplified GPR section. 
%
%       Author : Andreas Tzanis
%                Department of Geophysics, 
%                University of Athens
%                atzanis@geol.uoa.gr
%
%  Copyright (C) 2009, Andreas Tzanis. All rights reserved.
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
%
if isempty(d), 
    erh = errordlg('No data - the array is empty!','GAININVDECAY : ERROR');
     uiwait(erh);   clear erh; 
     dout = [];
     return
end

% Work arrays for fitting attenuation models, global within the scope of
% "getattenuation"
global tt decaycurve

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
ht = abs( ht(1:ns,:) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
attn       = [median( ht, 2 )  mean( ht, 2 ) ];
curve      = 1;
decaycurve = attn(:,curve);
tt    = tt2w;
warning('off','MATLAB:rankDeficientMatrix');
Tol   = 1.0e-8;
order = 3;
while 1==1
    % Get best fitting multi-exponential attenuation model
    fitoptions = optimset('Display','off','MaxFunEvals',2000*order, ...
        'MaxIter',2000*order, 'TolX',100*Tol,'TolFun',Tol');
    a0   = 0:1:order;
    a    = fminsearch(@expow,a0,fitoptions);
    ll = max(size(decaycurve));
    A = zeros(ll, size(a,2));
    for j = 1:size(a,2);
        A(:,j) = exp(-a(j)*tt2w(:));
    end;
    c = A\log10(decaycurve(:));
    expower = 10.^(A*c);

    % Display the results
    invdfig = findobj('tag','invdecayfigure');
    if ~isempty(invdfig),
        figure(invdfig);
        clf;
        figuretools;
        matgprwindows('updateadd');
    else
        invdfig = figure('name','DECAY MODEL',...
            'tag','invdecayfigure', ...
            'numbertitle','off', ...
            'menubar','none', ...
            'CreateFcn',[' figuretools; matgprwindows(''setup''); ' ...
            'matgprwindows(''updateadd'');'],...
            'DeleteFcn',['clear global tt decaycurve; ' ...
            'matgprwindows(''updateremove'', ' ...
            'findobj(''tag'',''invdecayfigure''));']);
    end

    semilogy(tt2w, attn(:,1), tt2w, attn(:,2), tt2w, expower, 'linewidth', 1 )
    xlabel('Time (ns)');  ylabel('Instantaneous Amplitude');
    grid;
    lines = get(gca,'children');
    lgd = legend(lines,...
        'Best fitting attenuation model', ...
        'Mean Attenuation', ...
        'Median Attenuation');
    set(lgd,'fontsize',10)

    if curve == 1, 
        altattn = 'Mean';  
    elseif curve == 2,
        altattn = 'Median';
    end
    nextthing = menu('SELECT OPTION','Proceed','Increase Order',...
        'Decrease Order',['Use ' altattn ' Attenuation'],'Cancel');

    switch nextthing
        case 1
            g    = 1./(expower/max(expower)) * ones(1, ntr);
            dout = d.*g;
            warning('on','MATLAB:rankDeficientMatrix');
            delete(invdfig);
            return
        case 2
            order = order + 1;
            continue
        case 3
            order = order - 1;
            continue
        case 4
            if curve == 1,
                curve      = 2;
            elseif curve ==2,
                curve      = 1;
             end
             decaycurve = attn(:,curve);
       case 5
            dout = [];
            warning('on','MATLAB:rankDeficientMatrix');
            delete(invdfig);
            return
    end
end

return
%
function er = expow(a)
global decaycurve tt
ll = max(size(decaycurve));
A = zeros(ll, size(a,2));
for j = 1:size(a,2);    
    A(:,j) = exp(-a(j)*tt(:)); 
end;
c = A\log10(decaycurve(:));
res = A*c;
% er = norm( log10(decaycurve(:)) - log10(res) );
er = norm( log10(decaycurve(:)) - res );
return