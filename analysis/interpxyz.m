function xyz = interpxyz(dx, ntr, markertr)
%
%  INTERPXYZ : Generate x, y and z trace coordinates with respect to 
%              a local (survey) reference system, by interpolating 
%              between the known coordinates of control (marker) traces. 
%              The X, Y and Z coordinates of the control traces are
%              passed as arguments or read from a disk (marker) file.
%              The traces of the input radargram must be equally spaced. 
%              The output x, y and z data is equally spaced in only one
%              of the x or y axes, (independent coordinate), chosen to be 
%              the longest monotonically varying dimension. However, the 
%              spacing of h = sqrt(x^2 + y^2) is roughly equal to the
%              trace spacing.
%
%      Usage : xyz = interpxyz(dx, ntr, markertr)
%
%     Inputs : 
%         dx : Trace spacing of the input radargram data. (empty if
%              data not equally spaced).
%        ntr : the number of traces in the input radargram.
%   markertr : [n x 4 ] matrix with the cooordinates of control (marker) 
%              traces: The columns respectively are the id numbers of 
%              control traces and their x, y and z coordinates in a
%              local frame of reference. 
%
%    Outputs :
%        xyz : Trace coordinates in a local frame of reference
%
%     Author : Andreas Tzanis, 
%              Department of Geophysics, 
%              University of Athens
%              atzanis@geol.uoa.gr
%              (C) 2005, Andreas Tzanis, all rights reserved
%

if nargin < 3,         % marker trace information not supplied
    erh = errordlg('Not enough input arguments','INTERPXYZ : ERROR');
    uiwait(erh);
    xyz = [];  
    return
end
if isempty(markertr),  % No marker traces
    erh = errordlg('Marker trace information does not exist! Aborting!',...
        'INTERPXYZ: ERROR');
    uiwait(erh);
    xyz = [];  
    return
end
if isempty(dx),        % Unequally spaced traces
    erh = errordlg('The traces are not equally spaced. Aborting!',...
        'INTERPXYZ: ERROR');
    uiwait(erh);
    xyz = [];  
    return
end
%%% Proceed with the analysis
sf = size(markertr); 
if sf(2)==4,  
    Markers   = markertr(:,1);
    lmtr      = length(Markers);
    Xmarker   = markertr(:,2);
    minXm     = min(Xmarker);
    Ymarker   = markertr(:,3);
    minYm     = min(Ymarker);
    Zmarker   = markertr(:,4);
% Initialize work variables
    X = [];     Y = [];      Z = [];
% Choose the longest monotonic coordinate X or Y to be the independent 
% coordinate. 
    xmono = 1;
    test=diff(Xmarker);
    ii = find(sign(test)>0);   jj = find(sign(test)<0);
    if ~isempty(ii) & ~isempty(jj)
        xmono = 0;
    end
    ymono = 1;
    test=diff(Ymarker);
    ii = find(sign(test)>0);   jj = find(sign(test)<0);
    if ~isempty(ii) & ~isempty(jj)
        ymono = 0;
    end
    if xmono & ymono,
        Hx = abs(max(Xmarker)-min(Xmarker));
        Hy = abs(max(Ymarker)-min(Ymarker));
        if Hx > Hy,
            indep_coord = 'X';
        else
            indep_coord = 'Y';
        end
    end
    if ~xmono & ymono,
        indep_coord = 'Y';
    end
    if xmono & ~ymono,
        indep_coord = 'X';
    end
    if ~xmono & ~ymono,
        answer = questdlg({'Neither x nor y are monotonic! ' ...
                'CONTINUE ONLY AT YOUR OWN RISK ... !!!'}, ...
                'INTERPXYZ: WARNING','Exit','Continue','Exit');
        if isequal(lower(answer),'exit'),
            xyz = [];
            return
        else
            Hx = abs(max(Xmarker)-min(Xmarker));
            Hy = abs(max(Ymarker)-min(Ymarker));
            if Hx > Hy,
                indep_coord = 'X';
            else
                indep_coord = 'Y';
            end
        end
    end
    disp(['INTERPXYZ > Independent coordinate is the ' indep_coord '-axis'])

% Independent coordinate is the x-axis
    if strcmp(indep_coord,'X'),       
% If traverse along the chosen indepenent coordinate runs opposite to 
% coordinate system axis, flip the marker data for the sake of 
% interpolation only!
        Xmarker       = Xmarker - minXm;
        Ymarker       = Ymarker - minYm;
        flipped = 0;
        if Xmarker(1) > Xmarker(length(Xmarker)),
            flipped = 1;
            Xmarker       = flipud(Xmarker)-minXm;
            Ymarker       = flipud(Ymarker)-minYm;
            Zmarker       = flipud(Zmarker);
        end
        delta_x = (max(Xmarker)-min(Xmarker))/(ntr-1);
        X       =  min(Xmarker) : delta_x : max(Xmarker);
        Y       = [];
% This computes a function y(x), equally spaced in x
        for i=1:ntr-1,
            ii=find(Xmarker <= X(i));
            j = ii(length(ii));
            slope = (Ymarker(j+1)-Ymarker(j))/(Xmarker(j+1)-Xmarker(j));
            intercept = Ymarker(j) - Xmarker(j)*slope;
            Y(i) = intercept + X(i)*slope;
        end
        slope = (Ymarker(lmtr)-Ymarker(lmtr-1))/...
            (Xmarker(lmtr)-Xmarker(lmtr-1));
        intercept = Ymarker(lmtr-1) - Xmarker(lmtr-1)*slope;
        Y(ntr) = intercept + X(ntr)*slope;
% This computes a function z(x), equally spaced in x
        Z = [];
        for i=1:ntr-1,
            ii=find(Xmarker <= X(i));
            j = ii(length(ii));
            slope = (Zmarker(j+1)-Zmarker(j))/(Xmarker(j+1)- Xmarker(j));
            intercept = Zmarker(j) - Xmarker(j)*slope;
            Z(i) = intercept + X(i)*slope;
        end
        slope = (Zmarker(lmtr)-Zmarker(lmtr-1))/...
            (Xmarker(lmtr)- Xmarker(lmtr-1));
        intercept = Zmarker(lmtr-1) - Xmarker(lmtr-1)*slope;
        Z(ntr) = intercept + X(ntr)*slope;
% If data was flipped, restore to original direction    
        if flipped, 
            X  = fliplr(X);
            Y  = fliplr(Y);
            Z  = fliplr(Z);
        end
        xyz = [X+minXm ;  Y+minYm ;  Z ]';
        disp('INTERPXYZ > Trace coordinate data created')
        return
    end

% Independent coordinate is the y-axis
    if strcmp(indep_coord,'Y'),
% If traverse along the chosen indepenent coordinate runs opposite to 
% coordinate system axis, flip the marker data to interpolate.
        Xmarker       = Xmarker - minXm;
        Ymarker       = Ymarker - minYm;
        flipped = 0;
        if Ymarker(1) > Ymarker(length(Ymarker)),
            flipped = 1;
            Xmarker       = flipud(Xmarker);
            Ymarker       = flipud(Ymarker);
            Zmarker       = flipud(Zmarker);
        end
% This computes a function x(y), equally spaced in y
        delta_y = (max(Ymarker)-min(Ymarker))/(ntr-1);
        Y       =  min(Ymarker) : delta_y : max(Ymarker);
        X       = [];
        for i=1:ntr-1
            ii=find(Ymarker <= Y(i));
            j = ii(length(ii));
            slope = (Xmarker(j+1)-Xmarker(j))/(Ymarker(j+1)-Ymarker(j));
            intercept = Xmarker(j) - Ymarker(j)*slope;
            X(i) = intercept + Y(i)*slope;
        end
        slope = (Xmarker(lmtr)-Xmarker(lmtr-1))/...
            (Ymarker(lmtr)-Ymarker(lmtr-1));
        intercept = Xmarker(lmtr-1) - Ymarker(lmtr-1)*slope;
        X(ntr) = intercept + Y(ntr)*slope;
% This computes a function z(y), equally spaced in y
        Z = [];
        for i=1:ntr-1,
            ii=find(Ymarker <= Y(i));
            j = ii(length(ii));
            slope = (Zmarker(j+1)-Zmarker(j))/(Ymarker(j+1)- Ymarker(j));
            intercept = Zmarker(j) - Ymarker(j)*slope;
            Z(i) = intercept + Y(i)*slope;
        end
        slope = (Zmarker(lmtr)-Zmarker(lmtr-1))/...
            (Ymarker(lmtr)- Ymarker(lmtr-1));
        intercept = Zmarker(lmtr-1) - Ymarker(lmtr-1)*slope;
        Z(ntr) = intercept + Y(ntr)*slope;
% If data was flipped, restore to original direction    
        if flipped, 
            X  = fliplr(X);
            Y  = fliplr(Y);
            Z  = fliplr(Z);
        end
        xyz = [X+minXm ;  Y+minYm ;  Z ]';
        disp('INTERPXYZ > Trace coordinate data created')
        return
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif sf(2)==1,
    erh = warndlg({'I can find only marker ID numbers but no x,y, z ' ...
        'data! Please, provide complete xyz information ' ...
        'for marker traces and try again!' }, ...
        ' INTERPXYZ: WARNING');
    uiwait(erh);
    xyz = [];  
    return
else
    erh = warndlg({'There are more than one columns in the variable '...
           '"markertr", but the structure is not standard! ' ...
           'Please edit and check!' }, ... 
           'INTERPXYZ: WARNING');
    uiwait(erh);
    xyz = [];  
    return
end
return
%%  The following code allows the generation of x,y and z data for raw
%%  radargrams recorded without survey wheels and such triggering
%%  facilities, in which traces are NOT equally spaced (as determined
%%  from the value of "dx"). In this case, the output x,y and
%%  z data are interpolated assuming constant recording rates bettween
%%  successive control points and is not equally spaced.
%%
%    if Scans_per_meter == 0,
%% First work out X, Y and Z between marked traces
%        for k = 1:lmtr-1,
%            dmark = Markers(k+1) - Markers(k);
%            xdmark = (Xmarker(k+1) - Xmarker(k))/dmark;
%            X = [X   Xmarker(k) : xdmark : Xmarker(k+1)-xdmark];
%            ydmark = (Ymarker(k+1) - Ymarker(k))/dmark;
%            Y = [Y   Ymarker(k) : ydmark : Ymarker(k+1)-ydmark];
%            zdmark = (Zmarker(k+1) - Zmarker(k))/dmark;
%            Z = [Z   Zmarker(k) : zdmark : Zmarker(k+1)-zdmark];
%        end
%        if xdmark==0,
%            X = Xmarker(1)+zeros(1,Markers(lmtr) - Markers(1) +1);
%        else 
%            X = [X X(length(X))+xdmark];
%        end
%        if ydmark==0,
%            Y = Ymarker(1)+zeros(1,Markers(lmtr) - Markers(1) +1);
%        else
%            Y = [Y Y(length(Y))+ydmark];
%        end
%        if zdmark==0,
%            Z = Zmarker(1)+zeros(1,Markers(lmtr) - Markers(1) +1);
%        else
%            Z = [Z Z(length(Z))+zdmark];
%        end
%% The next lines extrapolate beyond the first and last marked trace,
%% assuming that the recording rates were similar to those in the first and
%% last marked intervals 
%%        xdmark = (Xmarker(2) - Xmarker(1))/(Markers(2) - Markers(1));
%%        X = [ X(1)-fliplr(cumsum(ones(1,Markers(1)-1)*xdmark))  X  ];
%%        xdmark = (Xmarker(lmtr) - Xmarker(lmtr-1))/...
%%           (Markers(lmtr) - Markers(lmtr-1));
%%        X = [ X  X(length(X))+cumsum(ones(1,ntr-Markers(lmtr))*xdmark) ];
%%        ydmark = (Ymarker(2) - Ymarker(1))/(Markers(2) - Markers(1));
%%        Y = [ Y(1)-fliplr(cumsum(ones(1,Markers(1)-1)*ydmark))  Y  ];
%%        ydmark = (Ymarker(lmtr) - Ymarker(lmtr-1))/...
%%            (Markers(lmtr) - Markers(lmtr-1));
%%        Y = [ Y  Y(length(Y))+cumsum(ones(1,ntr-Markers(lmtr))*ydmark) ];
%%        Z = [ ones(1,Markers(1)-1)*Z(1)  Z ...
%%            ones(1,ntr-Markers(lmtr))*Z(length(Z))];
%%    
%        xyz = [ X' Y' Z'];
%        return          % unequally spaced profile - program returns here
%    end
