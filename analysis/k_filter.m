function  [dkf, Ws] = k_filter(tt, d, dx, filtertype, design_option)
%
%    K_FILTER : Design and apply FIR frequency filters. This routine is
%               exactly the same as F_FILTER, only that it operates
%               line-wise instead of column-wise, thus filtering
%               wavenumbers instead of traces (time series).
%
%       Usage : [dkf, Ws] = k_filter(tt,d,dx,filtertype,design_option) 
%
%      Inputs : 
%          tt : Time co-ordinates 
%           d : 2-D array of distance vs time input data ;   
%          dx : Trace spacing (constant along the scan line);
%  filtertype : flag determining the type of the filter;
%               =1  Low-pass filter
%               =2  High-pass filter
%               =3  Band-pass filter
%               =4  Band-stop filter
%design_option: keyword determining how the filter will be designed.
%               = 'TESTTR' : Use a test line. Point and click on some test
%                            line, then define the cutoff wavenumber(s)
%                            using the spectrum of the test line as a
%                            template, also by pointing and clicking. 
%               = 'MEANTR' : Compute the mean (average) scanline, then
%                            define the cutoff wavenumber(s) by pointing
%                            and clicking, using the spectrum of the mean
%                            line as template. 
%               = 'TYPEIN' : Enter predefined filter characteristics
%                            through dialog boxes.
%
%     Outputs : 
%        dkf  : 2-D array of filtered data
%         Ws  : The cutoff frequency (filtertype = 1 or 2), or, the pass /
%               stop bands (filtertype = 3 or 4).
%
%   Requires  :  fir_f1.m, fir_f2.m and checkcomma.m
%
%     Author  :  Andreas Tzanis,
%                Department of Geophysics, 
%                University of Athens
%                atzanis@geol.uoa.gr
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%

if filtertype < 1 | filtertype >4 | isstr(filtertype), 
    erh = errordlg('Wrong filter type. Aborting!','K_FILTER : ERROR');
    uiwait(erh);
    return
end;
%%% ns and ntr is Samples_per_scan and No_traces respectively
[ns, ntr] = size(d);
kn=1/(2*dx);
dk=kn/((ntr)/2);
k=dk:dk:kn; 

if strcmp(design_option,'TESTTR') || strcmp(design_option,'MEANTR'),
    
    if strcmp(design_option,'TESTTR'), 
%%%%%%%%   Get a scanline trace to design filter 
        datafig = findobj('tag','datafigure');
        if ~ishandle(datafig),
            erh = errordlg('No data figure for the filter design interface',...
                'K_FILTER : ERROR');
            uiwait(erh);
            dkf = [];
            return
        end
%%%% display the help message
        message     = msgbox('Please click on a test line',...
            'K_FILTER : HELP');
        messagepos  = get(message,'position');
        set(message,'pos',[20 60 150 messagepos(4)]);
        pause(0.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(datafig)
        pause(0.5)
        [xtr,ytr]=getpoint(' (m)', ' (ns)');        %[kcutoff,a]=ginput(1);
        xnum = round((ytr-tt(1))/(tt(2)-tt(1)));
        x  = d(xnum,:);
        x  = x-mean(x);
     
    elseif strcmp(design_option,'MEANTR'),
         x = mean(d);
         xnum = 0;
    end
    ftx=fft(x); 
%%%%%%%%  Design the filter on the test spectrum
    reply = 'No';
    while strcmp(reply,'Yes')~=1,
        if exist('trspfig') & ishandle(trspfig),
            figure(trspfig)
        else
            trspfig = figure('name','Select Cutoff Wavenumbers',...
                'numbertitle','off');
        end
        plot(k,abs(ftx(1:floor(ntr/2)))); set(gca,'nextplot','add');
        axis tight
        xlabel('Wavenumber (m^-1)');  
        ylabel(['Amplitude spectrum of line ' num2str(xnum) ]); 
        if filtertype == 1,
            [kcutoff,a]=getpoint(' m^(-1)',' ');    %[kcutoff,a]=ginput(1);
            kcutoff(1) = checkf(kcutoff(1),dk,kn);
            Ws=kcutoff(1)/kn;
            txt0 = ['Cutoff wavenumber = ' num2str(kcutoff,5) ' m^-^1'];
        elseif filtertype == 2,
            [kcutoff,a]=getpoint(' m^(-1)',' ');    %[kcutoff,a]=ginput(1);
            kcutoff(1) = checkf(kcutoff(1),dk,kn);
            Ws=kcutoff(1)/kn;
            txt0 = ['Cutoff wavenumber = ' num2str(kcutoff,5) ' m^-^1'];
        elseif filtertype == 3 || filtertype == 4,
            for i=1:2 
                [kcutoff,a]=getpoint(' m^(-1)',' ');%[kcutoff,a]=ginput(1);
                kcutoff(1) = checkf(kcutoff(1),dk,kn);
                Ws(i) = kcutoff(1)/kn;
                ylimit = get(gca,'ylim');
                plot([kcutoff(1) kcutoff(1)],[0 ylimit(2)], ':r', ...
                    'linewidth',2)
                if i==1,
                    text1 = ['k_1=' num2str(kcutoff(1),5) ];
                    text(kcutoff(1),0.9*ylimit(2),text1,'fontsize',9, ...
                        'fontwei','bold','horizontalal','right');
                end
                if i==2,
                    text2 = ['k_2=' num2str(kcutoff(1),5) ];
                    text(kcutoff(1),0.9*ylimit(2),text2,'fontsize',9, ...
                        'fontwei','bold','horizontalal','left');
                end
                drawnow
            end
           txt0 = ['Cutoff wavenumbers : ' text1 ' m^-^1 - ' text2 ' m^-^1'];
       end
% Compute the filter coefficients
       if filtertype == 1,
           b = fir_f1(floor(3*ntr/4), Ws);
       elseif filtertype == 2,
           b = fir_f1(floor(3*ntr/4), Ws, 'high');
       elseif filtertype ==3, 
           b = fir_f1(floor(3*ntr/4), Ws);
       elseif filtertype == 4,
           b = fir_f1(floor(3*ntr/4), Ws, 'stop' );
       end
% Test the filter
       if ~isempty(b),        %  An error in FIR_F1/2 produced an empty  b 
           [rb,cb]=size(b);
           if rb ~= 1, 
               b=b'; 
           end;
           b = [b zeros(1,ntr-length(b))];
           H     = fft(b); 
           ftxf  = ftx.*H;
           ftxf  = ftxf.*conj(H);
           figure(trspfig); 
           cla
           plot(k,abs(ftx(1:floor(ntr/2)))/max(abs(ftx)), ...
               k,abs(ftxf(1:floor(ntr/2)))/max(abs(ftx)), ...
               k,abs(H(1:floor(ntr/2))),'r');
           xlabel('Wavenumber (m^-^1)');  
           ylabel(['Amplitude spectrum of line ' num2str(xnum) ]); 
           title(txt0);
           axis tight;
           drawnow;
           reply = questdlg('Is this OK? ');
           if strcmp(reply,'Cancel')==1,
               delete(trspfig); 
               dkf=[];    
               return
           end
           if strcmp(reply,'No'),
               cla
           end
       else
           reply = 'No';
           cla
       end
   end
   delete(trspfig); 
   if exist('message') & ishandle(message),
       delete(message)
   end
end

%%% Cutoff wavenumbers are predefined and directly typed in 
if strcmp(design_option,'TYPEIN'),
    if filtertype == 1 | filtertype == 2, 
        answer = inputdlg('给出截止波数',...
            'K_FILTER: REQUEST',1);
        if isempty(answer), 
            dkf = []; 
            return, 
        end;
        lb = checkcomma(answer);        % decimals will be correct even if , is typed instead of .
        Ws = str2num(lb(1))/kn;
    end;
    if filtertype == 3 | filtertype == 4, 
        answer = inputdlg({'给出较低的截止波数 (m^-^1)',...
            '给出较高的截止波数 (m^-^1)'}, ...
            'K_FILTER: REQUEST',1);
        if isempty(answer), 
            dkf = []; 
            return, 
        end;
        lb = checkcomma(answer);        
        Ws = [str2num(lb(1,:)) str2num(lb(2,:))]/kn;
    end;
    if filtertype == 1,
        b = fir_f1(floor(3*ns/4), Ws);
    elseif filtertype == 2,
        b = fir_f1(floor(3*ns/4), Ws, 'high');
    elseif filtertype == 3,
        b = fir_f1(floor(3*ns/4), Ws);
    elseif filtertype == 4,
        b = fir_f1(floor(3*ns/4), Ws, 'stop' );
    end;
    if isempty(b),
        erh = errordlg('Error in filter computations! Please try again!',...
                'K_FILTER : ERROR');
        uiwait(erh);
        return
    end
    [rb,cb]=size(b);
    if rb ~= 1, 
        b=b'; 
    end;
    b = [b zeros(1,ntr-length(b))];
end

%%%% Filter in the frequency domain twice to preserve phase
hw = helpdlg('Processing... Please wait ...');
drawnow

H     = fft(b); 
H     = transpose(H) * ones(1, ns);   % filter already in k-domain - 
                                      % create ntr x ns matrix operator
ftx  = fft(d');                       % transform data to k-domain
ftx  = ftx.*H;                        % filter forwards
ftx  = ftx.*conj(H);                  % filter backwards
dkf   = real(ifft(ftx));              % invert to space domain
dkf   = dkf';

close(hw);
return

function fo = checkf(fi,dx,xmax)
if abs(fi) < abs(dx),
    fo = dx;
elseif abs(fi) > abs(xmax),
    fo = xmax;
else
    fo = fi;
end
