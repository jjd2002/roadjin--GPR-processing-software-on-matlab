function [header,A,samples,traces,samplingfrequency,DI] = loadhcd
%LOADRD3 loads radar profiles from a ramac format file. 
%
%Usage:%A = loadhcd;  opens an "open file" dialog box to interactively choose
%file. Matrix A is created. rows of A have common two-way-travel-time,
%columns of A have common profile station.  
%
%[A,samples,traces] = loadhcd;  will return matrix A and variables
%"sapmples" and "traces", which contain the dimensions of A.

%Copyright: (c) Sixin Liu, 2024
%Jilin Univ

%ask for file to load - automatically determines header file
[dat,pathname]= uigetfile ('*.hcd', 'Load Ramac format radar file');
name=double(dat);
chrctrs=size(name);
header=char([name(1:chrctrs(2)-2),97,100]);

%open header file - find number of samples and number of traces
fid=fopen([pathname,header],'r');
for k=1:50   % Read the Header
   line=fgetl(fid);   
   if ~isempty(findstr('SAMPLES:',line))      
      pt=findstr('SAMPLES:',line);
      samples=str2num(line(pt+8:length(line)));
   end
   if ~isempty(findstr('TRACE NUMBER:',line))
      pt=findstr('TRACE NUMBER:',line);
      traces=str2num(line(pt+13:length(line)));
   end
      if ~isempty(findstr('FREQUENCY:',line))
      pt=findstr('FREQUENCY:',line);
      samplingfrequency=str2num(line(pt+10:length(line)));
   end
      if ~isempty(findstr('USER DISTANCE INTERVAL:',line))
      pt=findstr('USER DISTANCE INTERVAL:',line);
      DI=str2num(line(pt+23:length(line)));
      end
      if ~isempty(findstr('DATA BIT:',line))
      pt=findstr('DATA BIT:',line);
      bit=str2num(line(pt+9:length(line)))
      end
end

fclose(fid);
%close header

%open data
fid=fopen([pathname,dat],'r');
if bit==16
A=fread(fid,[samples,traces],'int16');
elseif bit==32
   A=fread(fid,[samples,traces],'int32');
end
    
fclose(fid);



% function [DATA, HEADER] = loadhcd
% % LOADHCD : Imports HCD format georadar data.
% %
% %   Usage : DATA = loadhcd;     [DATA, HEADER] = loadhcd;
% %
% % RETURNS : A structure DATA containing the GPR data and all other necessary 
% %           parameters and, optionally, a structure HEADER containing the 
% %           header information of the imported HCD data file
% %
% % REQUIRES: MATLAB built-in functions (no external dependencies)
% %
% % Author : [Your Name]
% % Date   : [Date]
% 
% % Initialize the DATA structure
% DATA = struct();
% 
% % Ask user to select the HCD file
% [dat, pathname] = uigetfile('*.hcd', 'Load HCD format radar file');
% if isequal(dat, 0)
%     disp('User canceled the file selection');
%     return;
% end
% 
% % Generate the header file name (assuming it has the same name but with a different extension)
% headerFile = [pathname, strrep(dat, '.hcd', '.hdr')];
% 
% % Open and read the header file
% try
%     fid = fopen(headerFile, 'r');
%     if fid == -1
%         error('Unable to open header file: %s', headerFile);
%     end
% 
%     HEADER = struct(); % Initialize HEADER structure
% 
%     while ~feof(fid)
%         line = fgetl(fid);
%         if ischar(line)
%             if contains(line, 'SAMPLES:')
%                 HEADER.samples = str2double(extractAfter(line, 'SAMPLES:'));
%             elseif contains(line, 'TRACE NUMBER:')
%                 HEADER.traces = str2double(extractAfter(line, 'TRACE NUMBER:'));
%             elseif contains(line, 'FREQUENCY:')
%                 HEADER.samplingfrequency = str2double(extractAfter(line, 'FREQUENCY:'));
%             elseif contains(line, 'USER DISTANCE INTERVAL:')
%                 HEADER.DI = str2double(extractAfter(line, 'USER DISTANCE INTERVAL:'));
%             elseif contains(line, 'DATA BIT:')
%                 HEADER.bit = str2double(extractAfter(line, 'DATA BIT:'));
%             end
%         end
%     end
%     fclose(fid);
% catch ME
%     disp('Error reading header file');
%     disp(ME.message);
%     return;
% end
% 
% % Open and read the data file
% try
%     fid = fopen([pathname, dat], 'r');
%     if fid == -1
%         error('Unable to open data file: %s', [pathname, dat]);
%     end
% 
%     if HEADER.bit == 16
%         DATA.A = fread(fid, [HEADER.samples, HEADER.traces], 'int16');
%     elseif HEADER.bit == 32
%         DATA.A = fread(fid, [HEADER.samples, HEADER.traces], 'int32');
%     else
%         error('Unsupported data bit depth: %d', HEADER.bit);
%     end
%     fclose(fid);
%     
%     % Assign values to DATA structure
%     DATA.ns = HEADER.samples;
%     DATA.ntr = HEADER.traces;
%     DATA.dt = 1000 / HEADER.samplingfrequency;
%     DATA.sigpos = []; % Assuming no specific information for this field
%     DATA.dx = HEADER.DI;
%     DATA.Antenna = []; % Assuming no specific information for this field
%     DATA.TxRx = []; % Assuming no specific information for this field
%     DATA.comments = ''; % Assuming no specific comments to include
%     
%     % Two-way traveltime
%     DATA.tt2w = 0 : DATA.dt : DATA.dt * (DATA.ns - 1);
%     DATA.zlab = 'Traveltime (ns)';
%     
%     % Scan Line axis
%     if isempty(DATA.dx)
%         DATA.x = 1 : 1 : DATA.ntr;
%         DATA.xlab = 'Scan Axis (# Traces)';
%     else
%         DATA.x = 0 : DATA.dx : (DATA.ntr - 1) * DATA.dx;
%         DATA.xlab = 'Scan Axis (meters)';
%     end
%     
% catch ME
%     disp('Error reading data file');
%     disp(ME.message);
% end
% 
% end