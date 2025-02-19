function OPD = Run_batch_job(IPD)
%
%  RUN_BATCH_JOB : 
%  Executes a batch job using the priority and parameters supplied by the 
%  ENVAR.QUEUE structure (prepared with Setup_batch_job.m).
%
%  Usage : OPD = Setup_batch_job(IPD)
%
%  Input : IPD, the current Input Data structure of matGPR
%
% Output : OPD, the Output Data structure of matGPR
%
%  Author : Andreas Tzanis, 
%           Department of Geophysics, 
%           University of Athens, 
%           atzanis@geol.uoa.gr
%
% Copyright (C) 2010, Andreas Tzanis. All rights reserved.
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
global decaycurve tt

% Inotialize the Output structure
OPD = discardprocdata; 
% Trap common errors
if isempty(IPD.d), 
   erh = errordlg('请先导入数据!', 'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end
% XPD is a temporary data structure cum work array
XPD = IPD;
hist = 'BATCH JOB: ';
for istep=1:20,                                      % 20 is a dummy number

%%% Adjust Time Zero
    if istep == ENVAR.QUEUE.sigpos.step,
        OPD = XPD;
        if isempty(ENVAR.QUEUE.sigpos.new_zero),
            ENVAR.QUEUE.sigpos.new_zero = round(ENVAR.QUEUE.sigpos.new_t0/XPD.dt);
        end
        ns  = size(XPD.d,1);
        newns = ns - ENVAR.QUEUE.sigpos.new_zero;
        % if mod(newns,2) ~=0, 
        %     ENVAR.QUEUE.sigpos.new_zero = ENVAR.QUEUE.sigpos.new_zero - 1;
        % end
        OPD.d    = XPD.d(ENVAR.QUEUE.sigpos.new_zero+1:ns,:);
        OPD.tt2w = 0 : XPD.dt : (newns-1)*XPD.dt;
        OPD.ns   = size(OPD.d,1);
        if ~isempty(XPD.DZThdgain),
            OPD.DZThdgain    = XPD.DZThdgain(ENVAR.QUEUE.sigpos.new_zero+1:ns);
        end
        text = ['Moved Time-zero (' num2str(ENVAR.QUEUE.sigpos.new_zero*XPD.dt,4) 'ns)'];
        hist = [hist text '; '];
        disp(['STEP ' num2str(ENVAR.QUEUE.sigpos.step) ': ' text])
        XPD = OPD;
%%% Trim Time Window
    elseif istep == ENVAR.QUEUE.trimtw.step,
        OPD = XPD;
        if isempty(ENVAR.QUEUE.trimtw.cutoff),
            ENVAR.QUEUE.trimtw.cutoff = round(ENVAR.QUEUE.trimtw.new_tn/XPD.dt);
        end
        OPD.d    = XPD.d(1:ENVAR.QUEUE.trimtw.cutoff,:); 
        OPD.tt2w = XPD.tt2w(1:ENVAR.QUEUE.trimtw.cutoff);
        OPD.ns   = size(OPD.d,1);
        if ~isempty(XPD.DZThdgain),
            OPD.DZThdgain    = XPD.DZThdgain(1:ENVAR.QUEUE.trimtw.cutoff);
        end
        text = ['Trimmed Time Window (last ' ...
            num2str((XPD.ns-OPD.ns)*XPD.dt,5) ' ns)'];
        hist = [hist text '; '];
        disp(['STEP ' num2str(ENVAR.QUEUE.trimtw.step) ': ' text])
        XPD = OPD;
%%% Remove Global Background        
    elseif istep == ENVAR.QUEUE.rmbackgr.step,
        OPD = XPD;
        OPD.d = rmbackgr(XPD.d);
        text = 'Removed Global Background';
        hist = [hist text '; '];
        disp(['STEP ' num2str(ENVAR.QUEUE.rmbackgr.step) ': ' text])
        XPD = OPD;
%%% Remove DC        
    elseif istep == ENVAR.QUEUE.rmdc.step,
        OPD = XPD;
        OPD.d = removedc(XPD.d);
        text = 'Removed DC';
        hist = [hist text '; '];
        disp(['STEP ' num2str(ENVAR.QUEUE.rmdc.step) ': ' text])
        XPD = OPD;
%%% Dewow
    elseif istep == ENVAR.QUEUE.dewow.step,
        OPD = XPD;
        OPD.d = dewow(XPD.d);
        text = 'Dewow''ed';
        hist = [hist text '; '];
        disp(['STEP ' num2str(ENVAR.QUEUE.dewow.step) ': ' text])
        XPD = OPD;
%%% Remove DZT header gain
    elseif istep == ENVAR.QUEUE.dzthgain.step,
        if ~isempty(XPD.DZThdgain),
            OPD = XPD;
            OPD.d = gainrmdzthdr(XPD.d, XPD.DZThdgain);
            text = 'Removed DZT header gain';
            hist = [hist text '; '];
            disp(['STEP ' num2str(ENVAR.QUEUE.dzthgain.step) ': ' text])
            XPD = OPD;
        else
            disp(['STEP '  num2str(ENVAR.QUEUE.dzthgain.step) ...
                ': This is not a GSSI (DZT) data set! No action taken!']);
        end
%%% Standard AGC
    elseif istep == ENVAR.QUEUE.stdagc.step,
        OPD = XPD;
        iwagc = floor(ENVAR.QUEUE.stdagc.wagc/XPD.dt);      % agc window in samples
        if iwagc < 1 || iwagc > XPD.ns/2,
            disp(['STEP ' num2str(ENVAR.QUEUE.stdagc.step) ': AGC window is ' ...
                'too short or too long. Aborting...!']);
            OPD = discardprocdata;
            continue
        end
        OPD.d = gainagc(XPD.d, XPD.dt, ENVAR.QUEUE.stdagc.wagc);
        text = 'Applied AGC';
        hist = [hist text '; '];
        disp(['STEP ' num2str(ENVAR.QUEUE.stdagc.step) ': ' text])
        XPD = OPD;
%%% Gaussian-tapered AGC 
    elseif istep == ENVAR.QUEUE.gtagc.step, 
        OPD = XPD;
        iwagc = floor(ENVAR.QUEUE.gtagc.wagc/XPD.dt);      % agc window in samples
        if iwagc < 1 || iwagc > XPD.ns/2,
            disp(['STEP ' num2str(ENVAR.QUEUE.stdagc.step) ': AGC window is ' ...
                'too short or too long. Aborting...!']);
            OPD = discardprocdata;
            continue
        end
        OPD.d = gaingagc(XPD.d, XPD.dt, ENVAR.QUEUE.gtagc.EPS, ENVAR.QUEUE.gtagc.wagc);
        text = 'Applied Gaussian-tapered AGC';
        hist = [hist text '; '];
        disp(['STEP ' num2str(ENVAR.QUEUE.gtagc.step) ': ' text])
        XPD = OPD;
%%% Inverse Power Decay   
    elseif istep == ENVAR.QUEUE.invpdecay.step,
        OPD = XPD;
        [ns, ntr] = size(XPD.d);
        g     = (XPD.tt2w.^ENVAR.QUEUE.invpdecay.pow)' * ones(1,ntr);
        OPD.d = XPD.d.*g;
        % scale factor
        s2 = (max(max(OPD.d))-min(min(OPD.d)));
        s1 = (max(max(XPD.d))-min(min(XPD.d)));
        scale = s1/s2;
        OPD.d = scale*OPD.d;
        text = [ 'Applied gain g(t)= ' num2str(scale,4) '*t^' num2str(ENVAR.QUEUE.invpdecay.pow,3)];
        hist = [hist text '; '];
        disp(['STEP ' num2str(ENVAR.QUEUE.invpdecay.step) ': ' text])
        XPD = OPD;
%%% Inverse Amplitude Decay        
      elseif istep == ENVAR.QUEUE.invadecay.step,
        OPD = XPD;
        [ns, ntr] = size(XPD.d);
        % compute the analytic signal  %
        ns3 = 3*ns;                         % zero padd traces by 3-fold
        ftr = fft(XPD.d,ns3);
        h   = zeros(ns3,1);
        if ns3 >0 && 2*fix(ns3/2)==ns3      % even
            h([1 ns3/2+1]) = 1;
            h(2:ns3/2) = 2;
        elseif ns3 > 0                      % odd
            h(1) = 1;
            h(2:(ns3+1)/2) = 2;
        end
        ht = ifft( ftr .* (h * ones(1,ntr) ));
        ht = abs( ht(1:ns,:) );
        if ENVAR.QUEUE.invadecay.model == 1,
            decaycurve = median( ht, 2 );
        end
        if ENVAR.QUEUE.invadecay.model == 2,
            decaycurve = mean( ht, 2 );
        end
        tt    = XPD.tt2w;
        warning('off','MATLAB:rankDeficientMatrix');
        Tol   = 1.0e-8;
        fitoptions = optimset('Display','off','MaxFunEvals',2000*ENVAR.QUEUE.invadecay.order, ...
            'MaxIter',2000*ENVAR.QUEUE.invadecay.order, 'TolX',100*Tol,'TolFun',Tol');
        a0   = 0:1:ENVAR.QUEUE.invadecay.order;
        a    = fminsearch(@expow,a0,fitoptions);
        ll = max(size(decaycurve));
        A = zeros(ll, size(a,2));
        for j = 1:size(a,2);
            A(:,j) = exp(-a(j)*XPD.tt2w(:));
        end;
        c = A\log10(decaycurve(:));
        expower = 10.^(A*c);
        g    = 1./(expower/max(expower)) * ones(1, ntr);
        OPD.d = XPD.d.*g;
        warning('on','MATLAB:rankDeficientMatrix');
        % Processing History
        text = 'Applied Gain by Inv. Amplitude Decay';
        hist = [hist text '; '];
        disp(['STEP ' num2str(ENVAR.QUEUE.invadecay.step) ': ' text])
        XPD = OPD;      
%%% Resample Time Axis       
     elseif istep == ENVAR.QUEUE.rstime.step, 
        OPD = XPD;
        ns  = size(XPD.d,1);
        r      = ENVAR.QUEUE.rstime.new_ns/ns;         % change in sampling rate
        OPD.dt = XPD.dt/r;                     % new sampling rate
        OPD.ns = ENVAR.QUEUE.rstime.new_ns;
        N      = 15;                           % 1/2 order of sinc summation in resample1
        OPD.d = resample1(XPD.d, ENVAR.QUEUE.rstime.new_ns, ns, N);
        OPD.tt2w = 0:OPD.dt:(OPD.ns-1)*OPD.dt; % New traveltime vector
        if ~isempty(XPD.DZThdgain),            % If DZT data, update gain
            OPD.DZThdgain = interp1(XPD.tt2w,XPD.DZThdgain,OPD.tt2w);
        end
        text = ['Resampled Time Axis, ' num2str(XPD.ns) ...
            '->' num2str(OPD.ns) ' samples'];
        hist = [hist text '; '];
        disp(['STEP ' num2str(ENVAR.QUEUE.rstime.step) ': ' text])
        XPD = OPD;      
%%% Resample Scan Axis        
      elseif istep == ENVAR.QUEUE.rsscan.step, 
        OPD = XPD;
        ntr    = size(XPD.d,2);
        r      = ENVAR.QUEUE.rsscan.new_ntr/ntr;         % change in sampling rate
        OPD.dx = XPD.dx/r;                       % new sampling rate
        OPD.ntr= ENVAR.QUEUE.rsscan.new_ntr;
        OPD.x = XPD.x(1) + 0:OPD.dx:(OPD.ntr-1)*OPD.dx;
        N      = 15;             % 1/2 order of sinc summation in resample1
        % To comply with resample1 (working along the fast dimension), the
        % radargram must be transposed
        OPD.d = resample1(XPD.d', ENVAR.QUEUE.rsscan.new_ntr, ntr, N);
        OPD.d = OPD.d';                         % restore transposed output
        %%%%%
        % Now resample marker trace id's, albeit roughly
        if ~isempty(OPD.markertr),
            OPD.markertr(:,1) = fix(XPD.markertr(:,1)*OPD.ntr/XPD.ntr) ;
            OPD.markertr(1,1) = XPD.markertr(1,1);
        end;
        % Must resample trace coordinates!
        % Case 1. No marker trace info (e.g. Pulse Ekko data, or imported
        % from SU or SEGY file formats).
        if ~isempty(OPD.xyz.Tx) && isempty(OPD.markertr),
            % Resample Source Antenna coordinates
            dn = XPD.ntr/OPD.ntr;
            x1 = interp1([1:1:XPD.ntr]',XPD.xyz.Tx(:,1),[dn:dn:dn*OPD.ntr]');
            y1 = interp1([1:1:XPD.ntr]',XPD.xyz.Tx(:,2),[dn:dn:dn*OPD.ntr]');
            z1 = interp1([1:1:XPD.ntr]',XPD.xyz.Tx(:,3),[dn:dn:dn*OPD.ntr]');
            OPD.xyz.Tx = [x1(:) y1(:) z1(:)];
            % Resample Receiver Antenna coordinates
            x1 = interp1([1:1:XPD.ntr]',XPD.xyz.Rx(:,1),[dn:dn:dn*OPD.ntr]');
            y1 = interp1([1:1:XPD.ntr]',XPD.xyz.Rx(:,2),[dn:dn:dn*OPD.ntr]');
            z1 = interp1([1:1:XPD.ntr]',XPD.xyz.Rx(:,3),[dn:dn:dn*OPD.ntr]');
            OPD.xyz.Rx = [x1(:) y1(:) z1(:)];
        end;
        % Case 2. Monostatic data imported from radar file with marker
        % trace info! 
        if ~isempty(OPD.xyz.Tx) && ~isempty(OPD.markertr),
            OPD.xyz.Tx = interpxyz(OPD.dx,OPD.ntr,OPD.markertr);
            OPD.xyz.Rx = OPD.xyz.Tx;
        end;
        % Procesing History        
        text = ['Resampled Scan Axis, ' num2str(XPD.ntr) ...
            '->' num2str(OPD.ntr) ' traces'];
        hist = [hist text '; '];
        disp(['STEP ' num2str(ENVAR.QUEUE.rsscan.step) ': ' text])
        XPD = OPD;      
%%% Equalize Traces
    elseif istep == ENVAR.QUEUE.equalize.step, 
        if isempty(ENVAR.QUEUE.equalize.base),
            disp(['STEP ' num2str(ENVAR.QUEUE.equalize.step) ...
                ': Equalization base NOT specified! Equalization aborted!']);
            OPD = XPD; 
            continue
        end
        OPD = XPD;
        [ns, ntr] = size(XPD.d);
        OPD.d = zeros(ns,ntr);
        % Equalize to base value
        for ieq = 1:ntr;
            ss = sum(abs(XPD.d(:,ieq)));
            factor = ENVAR.QUEUE.equalize.base/ss;
            OPD.d(:,ieq) = XPD.d(:,ieq)*factor;
        end;
        text = 'Equalized traces';
        hist = [hist text '; '];
        disp(['STEP ' num2str(ENVAR.QUEUE.equalize.step) ': ' text])
        XPD = OPD;      
%%% FID Frequency Filtering        
       elseif istep == ENVAR.QUEUE.firff.step,
        OPD = XPD;
        [ns, ntr] = size(XPD.d);
        nf = floor(0.75*ns);               
        fn = 1/(2*XPD.dt*0.001);        % Nyquist in fake GHz!
        if ENVAR.QUEUE.firff.filtertype == 1,
            Ws = ENVAR.QUEUE.firff.f1/fn;
            b = fir_f1(nf, Ws);
        elseif ENVAR.QUEUE.firff.filtertype == 2,
            Ws = ENVAR.QUEUE.firff.f1/fn;
            b = fir_f1(nf, Ws, 'high');
        elseif ENVAR.QUEUE.firff.filtertype == 3,
            Ws = [ENVAR.QUEUE.firff.f1 ENVAR.QUEUE.firff.f2]/fn;
            b = fir_f1(nf, Ws);
        elseif ENVAR.QUEUE.firff.filtertype == 4,
            Ws = [ENVAR.QUEUE.firff.f1 ENVAR.QUEUE.firff.f2]/fn;
            b = fir_f1(nf, Ws, 'stop' );
        end;
        if isempty(b),
            disp(['STEP ' num2str(ENVAR.QUEUE.firff.step) ': Error in '...
                'Frequency Filter Computations! Aborting...']);
            continue
        end
        b = [b; zeros(ns-length(b),1); zeros(ns,1)]; %#ok<AGROW>
        %%%% Filter in the frequency domain twice to preserve phase
        hw = helpdlg('Processing... Please wait ...');
        H     = fft(b);                         % Filter in frequency domain
        H     = H * ones(1, ntr);               % Cast filter into a matrix operator
        fttr  = fft([XPD.d; zeros(ns,ntr)]);    % transform data
        fttr  = fttr.*H;                        % filter forwards
        fttr  = fttr.*conj(H);                  % filter backwards
        OPD.d   = real(ifft(fttr));             % invert to time domain
        OPD.d   = OPD.d(1:ns,:);
        close(hw);
        clear H fttr
        if ENVAR.QUEUE.firff.filtertype == 1,
            text = ['Low-Pass Filtered at ' ...
                num2str(Ws(1)*fn) ' MHz'];
        elseif ENVAR.QUEUE.firff.filtertype == 2,
            text = ['High-Pass Filtered at ' ...
                num2str(Ws(1)*fn) ' MHz'];
        elseif ENVAR.QUEUE.firff.filtertype == 3,
            text = ['Band-Pass ltered, ['...
                num2str(Ws(1)*fn,4) '-' num2str(Ws(2)*fn,4) '] MHz'];
        elseif ENVAR.QUEUE.firff.filtertype == 4,
            text = ['Band-Stop Filtered, ['...
                num2str(Ws(1)*fn,4) '-' num2str(Ws(2)*fn,4) '] MHz'];
        end
        hist = [hist text '; '];
        disp(['STEP ' num2str(ENVAR.QUEUE.firff.step) ': ' text])
        XPD = OPD;
        
    end
end

%%% Case when the batch job queue is empty (not set up)!
if isempty(OPD.d),
    disp('RUN BATCH JOB > No O/P data returned! Please check the job queue!');
    return
end


%%% Update processing history
OPD = XPD;
iss = size(OPD.history,1);
OPD.history(iss+1,1) = cellstr(hist);

% Display O/P data
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab);

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

