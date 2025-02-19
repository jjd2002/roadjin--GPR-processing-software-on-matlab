function [di,ntri,dx,markertr,scanline] = ...
    markerinterp(d, markertr, Inpname, Infname)
%
%    MARKERINTERP : Marker Interpolation from GPR traces equally spaced in 
%                   time to traces equally spaced in offset. 
%                 * If the section was recorded in a direction reverse with
%                   respect to the positive axes of coordinate system, (for
%                   instance when the survey is conducted in meandering
%                   forward - backward sense), the section is flipped, as
%                   also are marker trace data. This streamlines the
%                   transformed section with the axes of the reference
%                   frame.  
%
%           Usage : [di,ntri,dx,markertr,scanline] = 
%                         markerinterp(d, markertr, Inpname, Infname)
%
%          Inputs : 
%               d : the 2-D GPR section to be transformed from equal time 
%                   spacing to equal-distance spacing
%        markertr : [n x 4 ] array with the cooordinates of control (marker) 
%                   traces: The columns respectively are the id numbers of 
%                   control traces and their x, y and z coordinates in a
%                   local frame of reference. 
%  Inpname,Infname: The path and name of the GPR data file
%
%         Outputs :
%              di : The transformed (equidistant trace) GPR section
%            ntri : The number of traces in the transformed radargram
%              dx : Trace spacing in the transformed radargram (in m)
%        markertr : [n x 4 ] array of updated marker trace information 
%        scanline : ntri row vector of equally spaced trace coordinates 
%                   along the Scan Line
%
%        Requires : checkcomma.m 
%
%          Author : Andreas Tzanis, 
%                   Department of Geophysics, 
%                   University of Athens
%                   atzanis@geol.uoa.gr
%                  (C) 2005, Andreas Tzanis, all rights reserved
%

%%% Inquire where the x,y,z data of marker traces (control points) 
%%% will come from.
ask = questdlg('Import marker trace data from disk file ?', ...
    'MARKERINTERP: REQUEST', 'No');
if isempty(ask) || strcmpi(ask,'cancel'),
    di = [];   ntri = [];   dx = [];   scanline = [];
    return 
end;
%%% If the xyz data should come from a marker file, open it and load 
if strcmpi(ask,'yes')==1,
    Markfile = [Infname(1:findstr(Infname,'.')) 'mrk'];
    [Markfile, Inpname]= uigetfile('*.mrk; *.MRK',...
        'Give Marker File name', [Inpname Markfile]); 
    if Markfile == 0,
        di = [];   ntri = [];   dx = [];   scanline = [];
        return
    end
    fid = fopen([Inpname Markfile],'r');
    markertr = [];
    while ~feof(fid), 
       textline = fgetl(fid);
       if isstr(textline), 
           markertr = [markertr; str2num(textline)]; 
        end;
    end
    fclose(fid);   
end
%%% Now proceed
sf = size(markertr); 
if sf(2)==4,  
    Markers  = markertr(:,1);
    Xmarker  = markertr(:,2) - min(markertr(:,2));
    Ymarker  = markertr(:,3) - min(markertr(:,3));
    Zmarker  = markertr(:,4);
%%% If traverse runs opposite to axes of coord system, flip the section as
%%% well as marker trace data => streamline the transformed section with
%%% the coordinate system. A must when data recorded ina a meandering
%%% forward - backward sense. 
    Hmarker = sqrt(Xmarker.^2 + Ymarker.^2);  % distance along scanline
    if Hmarker(1) > Hmarker(length(Hmarker)), % section run in reverse
        d = fliplr(d);                        % flip section 
        ntr = size(d,2);
        Markers  = flipud((ntr+1) - Markers ); % flip marker information
        Hmarker  = flipud(Hmarker);
        Xmarker  = flipud(Xmarker);
        Ymarker  = flipud(Ymarker);
        Zmarker  = flipud(Zmarker);
    end
%%% Construct vectors with approximate, unequal spacings between 
%%% the traces of the input data. Still, we need to assume that the 
%%% recording rate is constant between marker intervals. 
    scanla = [];    Xa = [];     Ya = [];     Za = [];      %work arrays
    lmtr  = length(Markers);
    for k = 1:lmtr-1
        dmark = abs(Markers(k+1) - Markers(k));
        xdmark = (Xmarker(k+1) - Xmarker(k))/dmark;
        Xa = [Xa   Xmarker(k) : xdmark : Xmarker(k+1)-xdmark];
        ydmark = (Ymarker(k+1) - Ymarker(k))/dmark;
        Ya = [Ya   Ymarker(k) : ydmark : Ymarker(k+1)-ydmark];
    end
    if xdmark==0,
        Xa = Xmarker(1)+zeros(1,Markers(lmtr) - Markers(1) +1);
    else 
        Xa = [Xa Xa(length(Xa)) + xdmark];
    end
    if ydmark==0,
        Ya = Ymarker(1)+zeros(1,Markers(lmtr) - Markers(1) +1);
    else
        Ya = [Ya (Ya(length(Ya)) + ydmark)];
    end
    scanla = sqrt(Xa.^2 + Ya.^2);
%%% Construct equaly spaced H-array for interpolation %%%%%%%%%%%%%%%%%%%%%
    ask = inputdlg('Give desired trace spacing in meters',...
        'MARKERINTERP : REQUEST',1);
    if isempty(ask),
        di = [];    ntri = [];   dx = [];   scanline = [];
        return
    end
    lb   = checkcomma(ask);
    dx   = str2num(lb(1,:));
    scanline  = min(Hmarker) : dx : max(Hmarker);
    ntri = length(scanline);
%%% For convenience in subsequent processing, make sure that the
%%% transformed data contain only even number of traces. 
    if rem(ntri,2) ~= 0,
        scanline = [scanline(1:ntri) scanline(ntri)+dx];
        ntri = ntri + 1;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif sf(2)==1,
    erh = warndlg({'I can find only marker ID numbers but no x,y, z data! ' ...
        'Please, provide complete marker trace information and try again!' }, ...
        ' MARKERINTERP : WARNING');
    uiwait(erh);
    di = [];   ntri = [];   dx = [];   scanline = [];
    return
else
    erh = warndlg({'There are more than one columns in the variable markertr, ' ...
           'but the structure is not standard! Please edit and check!' }, ... 
           ' MARKERINTERP: WARNING');
    uiwait(erh);
    di = [];   ntri = [];   dx = [];   scanline = [];
    return
end
%%% Interpolate to equal trace spacing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hw = waitbar(0,'Marker interpolation in progress, Please wait');
ns = size(d,1);
di = zeros(ns, ntri);
for i=1:ns
    dummy  = d(i,Markers(1):Markers(lmtr));
    dummyi = interp1(scanla,dummy,scanline,'pchip');
    di(i,:) = dummyi;
    waitbar(i/ns,hw);
end
close(hw)
%%% Update id's of Marker Traces
newMarkers = [];
for i=1:size(Markers,1); 
    ii=find(abs(scanline-Hmarker(i)) <= dx); 
    newMarkers = [newMarkers; ii(length(ii))]; 
end
markertr = [newMarkers Xmarker+min(markertr(:,2)) ...
    Ymarker+min(markertr(:,3)) Zmarker];
return
% END FUNCTION markerinterp