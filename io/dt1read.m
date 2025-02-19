function [data,trheaders] = dt1read(filename)
%
% DT1READ: Read bistatic PULSE EKKO GPR data for binary file .DT1 (ignore
% the header file .HD)
%
%  Usage : [data, trheaders] = dt1read(filename);
%
%   Input : 
%filename : The path and name of the DT1 data file
%
%   Output:
%   data  : The imported 2-D GPR section 
%trheaders: Array of data structures with the trace headers of the GPR data
%
% Requires: readdt1traceheader.m, dt1struct.m, dimstruct.m
%
% Author: 
% Luca Baradello,
% lbaradello@ogs.trieste.it
%
% Acknowledgements: 
% Lian Zhao for info
%
 

fid = fopen(filename,'rb');    
dt1 = dt1struct;     
h   = dimstruct;            
fseek(fid,h.samples,-1); 
samples = fread(fid,1,dt1.samples);         
tracedim= samples*2+128  ;                                  
fseek(fid,-tracedim,1);
max_traces = fread(fid,1,dt1.traces);
i = 1;
fseek(fid,0,'bof');
for j = 1:max_traces,
     i = i + 1;  
     trheaders(i-1) = readdt1traceheader(fid,j); 
     trace = fread(fid,samples,dt1.trace);
     data(:,i-1) = trace(:); 
     position = tracedim*j;
     fseek(fid,position,-1);
end 
fclose(fid);

function dt1 = dt1struct;
dt1.traces= 'float';
dt1.position= 'float';
dt1.samples= 'float';
dt1.topo= 'float';
dt1.x1= 'float';
dt1.bytes= 'float';
dt1.trac_num= 'float';
dt1.stack= 'float';
dt1.window= 'float';
dt1.x2= 'float';
dt1.x3= 'float';
dt1.x4= 'float';
dt1.x5= 'float';
dt1.x6= 'float';
dt1.x_rec= 'float';
dt1.y_rec= 'float';
dt1.z_rec= 'float';
dt1.x_tra= 'float';
dt1.y_tra= 'float';
dt1.z_tra= 'float';
dt1.time_zero= 'float';
dt1.zero= 'float';
dt1.x7= 'float';
dt1.time= 'float';
dt1.x8= 'float';
dt1.com0= 'float';
dt1.com1= 'float';
dt1.com2= 'float';
dt1.com3= 'float';
dt1.com4= 'float';
dt1.com5= 'float';
dt1.com6= 'float';
dt1.com7= 'float';
dt1.trace= 'short';

function [dim]=dimstruct;
dim.tracl= 0;
dim.traces= 0;
dim.position= 4;
dim.samples= 8;
dim.topo= 12;
dim.x1= 16;
dim.bytes= 20;
dim.trac_num= 24;
dim.stack= 28;
dim.window= 32;
dim.x2= 36;
dim.x3= 40;
dim.x4= 44;
dim.x5= 48;
dim.x6= 52;
dim.x_rec= 56;
dim.y_rec= 60;
dim.z_rec= 64;
dim.x_tra= 68;
dim.y_tra= 72;
dim.z_tra= 76;
dim.time_zero= 80;
dim.zero= 84;
dim.x7= 88;
dim.time= 92;
dim.x8= 96;
dim.com0= 100;
dim.com1= 104;
dim.com2= 108;
dim.com3= 112;
dim.com4= 116;
dim.com5= 120;
dim.com6= 124;
dim.com7= 128;
      
function trheader = readdt1traceheader(fid,j)
trheader.traces= fread(fid,1,'float');
trheader.position= fread(fid,1,'float');
trheader.samples= fread(fid,1,'float');
trheader.topo= fread(fid,1,'float');
trheader.x1= fread(fid,1,'float');
trheader.bytes= fread(fid,1,'float');
trheader.trac_num= fread(fid,1,'float');
trheader.stack= fread(fid,1,'float');
trheader.window= fread(fid,1,'float');
trheader.x2= fread(fid,1,'float');
trheader.x3= fread(fid,1,'float');
trheader.x4= fread(fid,1,'float');
trheader.x5= fread(fid,1,'float');
trheader.x6= fread(fid,1,'float');
trheader.x_rec= fread(fid,1,'float');
trheader.y_rec= fread(fid,1,'float');
trheader.z_rec= fread(fid,1,'float');
trheader.x_tra= fread(fid,1,'float');
trheader.y_tra= fread(fid,1,'float');
trheader.z_tra= fread(fid,1,'float');
trheader.time_zero= fread(fid,1,'float');
trheader.zero= fread(fid,1,'float');
trheader.x7= fread(fid,1,'float');
trheader.time= fread(fid,1,'float');
trheader.x8= fread(fid,1,'float');
trheader.com0= fread(fid,1,'float');
trheader.com1= fread(fid,1,'float');
trheader.com2= fread(fid,1,'float');
trheader.com3= fread(fid,1,'float');
trheader.com4= fread(fid,1,'float');
trheader.com5= fread(fid,1,'float');
trheader.com6= fread(fid,1,'float');
      
