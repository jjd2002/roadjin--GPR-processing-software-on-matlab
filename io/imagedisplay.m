function imagedisplay(Data,x,t,xlab,tlab)
%
% IMAGEDISPLAY  : plots image of data in d"
%
%  Usage        :  imagedisplay(d,x,t,xlab,ylab)
%
%        INPUTS :
%          Data : the data, [ns , ntr]
%             x : X axis, [1:ntr]
%             t : Time axis, [1:ns]
%    xlab, tlab : x- and time axis labels respectively.
%
%      REQUIRES : imagecolors.m
%
%        Author : Andreas Tzanis, 
%                 Department of Geophysics, 
%                 University of Athens
%       Created : 23 August 2003

[ns,ntr] = size(Data);
if nargin ==1,
    t=[1:1:ns];
    x=[1:1:ntr];
end
if nargin==2,
    y=[1:1:ns];
end

imagesc(x, t, Data);

if nargin == 5 ,
    xlabel(xlab);
    ylabel(tlab);
end
if nargin >= 4,
    xlabel(xlab);
end

colorbar('h');

return        
