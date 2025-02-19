function markerdisplay(d,x,markertr)
%
% MARKERDISPLAY : Points the location of Marker Traces on the "GPR DATA"
%                 figure 
%
%        Usage  : markerdisplay(d, x, markertr)
%
%            d  : The GPR section of size [ns, ntr] 
%            x  : [1 x ntr], vector of the horizontal dimension in data
%                 units, usually distance.  
%      markertr : The locations of the Marked Traces
%
%        Author : Andreas Tzanis
%                 Department of Geophysics, 
%                 University of Athens
%                 atzanis@geol.uoa.gr
%                 (C) 2005, Andreas Tzanis, all rights reserved.
%
% 
if isempty(markertr),
    erh = errordlg('There is no marker information to display!', ...
        'MARKERDISPLAY : ERROR');
    uiwait(erh);
    return
end
[ns,ntr]=size(d);
ii = find(markertr > ntr);
if ~isempty(ii),
    erh = errordlg(['The number of traces has changed. ' ... 
        'Cannot display the original marker info!'], ...
        'MARKERDISPLAY : ERROR');
    uiwait(erh);
    return
end
% Get data figure handle
datafig = findobj('tag','datafigure');
if ~isempty(datafig),
    figure(datafig);
end
set(gca,'nextplot','add');
plot(x(markertr(:,1)),zeros(length(markertr),1),'kv',...
    'markersize',5,'markerfacecolor','k');
return

    