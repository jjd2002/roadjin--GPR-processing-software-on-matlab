function dg = gainrmdzthdr( d, DZThdgain )
%
% GAINRMDZTHDR : Remove gain from GSSI (DZT) data using the range-gain
%                information supplied in the header.  
%
%        Usage : dg  = rmdzthdrgain(d, DZThdgain )
%   
%       Inputs : d           is the GPR section
%                DZThdgain   is the DZT header gain function in db 
%
%       Output : dg          is the output section (same size as d);
%
%       Author : Andreas Tzanis, 
%                Department of Geophysics, 
%                University of Athens
%                atzanis@geol.uoa.gr
%                (C) 2005, Andreas Tzanis, all rights reserved
%

[ns,ntr] = size(d); 
g  = 10.^(DZThdgain/20);    % Expand header gain from db to real amplitude
g  = g' * ones(1,ntr);          
dg = d./g ;                 % Remove header gain

return
