% Script for estimating modal decomposition of velocity from moorings.
% Vertical mode shapes are based on OCCA v.2 salinity and temperature.
%
% C. Wortham

clear all
close all

% uses BTTmodes.m function here:
path(path,'~worthamc/Documents/research/scripts/VertStruct')

% OCCA climatology is here:
OCCAdir = '/Users/stallion/Documents/research/data/OCCA_1x1_v2/annual/';

fno = 0;

%% load mooring data: lotus 1
% for lotus1 (May 1982 to April 1983)
moorname = 'lotus 1';

% lotus mooring is here
MOORdir = '~worthamc/Documents/research/data/currentmetersOTHER/lotus/';

lat = 33.9;
lon = -70.0 + 360;
dt = 1/96; % days
starttime = datenum(1982,5,11,12,07,30);
endtime = datenum(1983,4,12,7,07,30);
depth_MOOR = -[127 129 178 228 278 328 427 775 1024 1521 4007];
badval = -999.9990;

% load data
load([MOORdir 'lotus1.mat']);
instruments = fieldnames(lotus1);
time = lotus1.(instruments{1});
u = nan(length(depth_MOOR),length(time));
v = nan(length(depth_MOOR),length(time));
for ii=1:length(depth_MOOR)
    u(ii,:) = lotus1.(instruments{2*ii});
    v(ii,:) = lotus1.(instruments{2*ii+1});
end
time = linspace(starttime,endtime,length(time));

% interpolate short gaps in velocity data
for ii = 1:length(depth_MOOR)
    missing = u(ii,:) == badval;
    longgap = strfind(missing,ones(1,30)); % looks for long data gaps
    if isempty(longgap)
        u_temp = u(ii,~missing);
        v_temp = v(ii,~missing);
        time_temp = time(~missing);
        u(ii,:) = interp1(time_temp,u_temp,time,'linear','extrap')';
        v(ii,:) = interp1(time_temp,v_temp,time,'linear','extrap')';
        disp([num2str(sum(missing)) ' data points interpolated'])
    end
end

% remove time mean velocity
ubar = mean(u,2); ubar = repmat(ubar,[1,length(time)]);
vbar = mean(v,2); vbar = repmat(vbar,[1,length(time)]);
u = u - ubar;
v = v - vbar;

% make units m/s
u = u/100;
v = v/100;

% drop instrument 7 (failed half way through)
depth_MOOR = depth_MOOR([1 2 3 4 5 6 8 9 10 11]);
u = u([1 2 3 4 5 6 8 9 10 11],:);
v = v([1 2 3 4 5 6 8 9 10 11],:);

% lowpass filter
% fN = 1/2/dt; % Nyquist frequency (cph)
% f1 = 1/24; % cutoff frequency
% b = fir1(20,f1/fN); % lowpass Hamming window filter (default)
% ufilt = nan(size(u));
% vfilt = nan(size(u));
% for ii = 1:length(depth_MOOR)
%     ufilt(ii,:) = filtfilt(b,1,u(ii,:));
%     vfilt(ii,:) = filtfilt(b,1,v(ii,:));
% end
% filter and decimate to twice daily
r = round(.5/dt);
dt = r*dt;
for ii = 1:length(depth_MOOR)
    ufilt(ii,:) = decimate(u(ii,:),r);
    vfilt(ii,:) = decimate(v(ii,:),r);
end

u = ufilt;
v = vfilt;

%% load mooring data: lotus 2
% for lotus1 (April 1983 to May 1984)
moorname = 'lotus 2';

% lotus mooring is here
MOORdir = '~worthamc/Documents/research/data/currentmetersOTHER/lotus/';

lat = 33.9;
lon = -70.0 + 360;
dt = 1/96; % days
starttime = datenum(1983,4,14,12,07,30);
endtime = datenum(1983,5,06,20,07,30);
depth_MOOR = -[98 198 248 298 348 398 448 498 598 748 998 2498 3998];
badval = -999.9990;

% load data
load([MOORdir 'lotus2.mat']);
instruments = fieldnames(lotus2);
time = lotus2.(instruments{1});
u = nan(length(depth_MOOR),length(time));
v = nan(length(depth_MOOR),length(time));
for ii=1:length(depth_MOOR)
    u(ii,:) = lotus2.(instruments{2*ii});
    v(ii,:) = lotus2.(instruments{2*ii+1});
end
time = linspace(starttime,endtime,length(time));

% drop end of file (all instruments stop)
stopind = 36750; % instruments start crapping out after this
u = u(:,1:stopind);
v = v(:,1:stopind);
time = time(1:stopind);

% interpolate short gaps in velocity data
for ii = 1:length(depth_MOOR)
    missing = u(ii,:) == badval;
    longgap = strfind(missing,ones(1,30)); % looks for long data gaps
    if isempty(longgap)
        u_temp = u(ii,~missing);
        v_temp = v(ii,~missing);
        time_temp = time(~missing);
        u(ii,:) = interp1(time_temp,u_temp,time,'linear','extrap')';
        v(ii,:) = interp1(time_temp,v_temp,time,'linear','extrap')';
        disp(['instrument ' num2str(ii) ', ' num2str(sum(missing)) ' data points interpolated'])
    end
end

% remove time mean velocity
ubar = mean(u,2); ubar = repmat(ubar,[1,length(time)]);
vbar = mean(v,2); vbar = repmat(vbar,[1,length(time)]);
u = u - ubar;
v = v - vbar;

% make units m/s
u = u/100;
v = v/100;

% drop instruments 1,2,3 (failed half way through)
depth_MOOR = depth_MOOR([4 5 6 7 8 9 10 11 12 13]);
u = u([4 5 6 7 8 9 10 11 12 13],:);
v = v([4 5 6 7 8 9 10 11 12 13],:);

% lowpass filter
% fN = 1/2/dt; % Nyquist frequency (cph)
% f1 = 1/24; % cutoff frequency
% b = fir1(50,f1/fN); % lowpass Hamming window filter (default)
% ufilt = nan(size(u));
% vfilt = nan(size(u));
% for ii = 1:length(depth_MOOR)
%     ufilt(ii,:) = filtfilt(b,1,u(ii,:));
%     vfilt(ii,:) = filtfilt(b,1,v(ii,:));
% end
% filter and decimate to twice daily
r = round(.5/dt);
dt = r*dt;
for ii = 1:length(depth_MOOR)
    ufilt(ii,:) = decimate(u(ii,:),r);
    vfilt(ii,:) = decimate(v(ii,:),r);
end

u = ufilt;
v = vfilt;

%% CRAP load mooring data: WOCE ACM12 VE/338
% (12 Jan 91 - 05 Dec 92)
% poor choice... it's in the Vema Channel, just south of the brazil basin.
% Flow is bottom intensified here, and OCCA hydrography doesn't cover the
% bottom of the mooring.

% WOCE ACM12 mooring data is here
MOORdir = '~worthamc/Documents/research/data/currentmetersOSU/Atlantic/WOCE/acm12/';

filenames = {'rcm01762.str' 'rcm01763.str' 'rcm01764.str' 'rcm01765.str' 'rcm01766.str' 'rcm01767.str' 'rcm01768.str'};

% get the data
for ii = 1:length(filenames);
    fid = fopen([MOORdir filenames{ii}],'r');
    % collect information from header
    headlines = textscan(fid,'%d header lines',1);headlines = headlines{:};
    datalines = textscan(fid,'%d data lines',1);datalines = datalines{:};
    textscan(fid,'%s',1); % skipping a line
    expname = textscan(fid,'Experiment name: %[^\n]',1);expname = char(deblank(expname{1}));
    moorname = textscan(fid,'Mooring name: %[^\n]',1);moorname = char(deblank(moorname{1}));
    latlong = textscan(fid,'Mooring position: %f deg %c, %f deg %c',1);
    if (latlong{2}=='N'), lat = latlong{1}; else lat = -latlong{1}; end
    if (latlong{4}=='E'), lon = latlong{3}; else lon = 360-latlong{3}; end
    instdepth = textscan(fid,'Instrument depth: %d meters',1);instdepth = instdepth{1};
    waterdepth = textscan(fid,'Seafloor depth: %d meters',1);waterdepth = waterdepth{1};
    instrument = textscan(fid,'Instrument type: %[^\n]',1);instrument = char(deblank(instrument{1}));
    accession = textscan(fid,'CMDB accession number: %d',1);accession = accession{1};
    for i=1:5
        textscan(fid,'%s',1); % skipping some lines
    end
    velunits = textscan(fid,'speed (%[^)]',1);velunits = char(velunits{1});
    % load the data
    fseek(fid,0,'bof'); %reset file position indicator to beginning
    if headlines==20 % setting number of fields based on header length
        fields='%f %f %f %f %f %f %f %f %f';
    elseif headlines==21
        fields='%f %f %f %f %f %f %f %f %f %f';
    elseif headlines==22
        fields='%f %f %f %f %f %f %f %f %f %f %f';
    elseif headlines==23
        fields='%f %f %f %f %f %f %f %f %f %f %f %f';
    else
        disp('check number of fields')
    end
    theData = textscan(fid, fields, datalines, 'headerlines', headlines);
    fclose(fid);
    
    % add to final data vectors
    depth_MOOR(ii) = -double(instdepth);
    stopind = 8325;
    u(ii,:) = theData{7}(1:stopind);
    v(ii,:) = theData{8}(1:stopind);
    hr = theData{1}(1:stopind);
    mm = mod(hr,100);
    hr = (hr-mm)/100;
    time = datenum(theData{4}(1:stopind)+1900,theData{3}(1:stopind),theData{2}(1:stopind),hr,mm,0);
end

dt = (time(2) - time(1)); % days

% remove time mean velocity
ubar = mean(u,2); ubar = repmat(ubar,[1,length(time)]);
vbar = mean(v,2); vbar = repmat(vbar,[1,length(time)]);
u = u - ubar;
v = v - vbar;

% make units m/s
u = u/100;
v = v/100;

% lowpass filter
% fN = 1/2/dt; % Nyquist frequency (cph)
% f1 = 1/12; % cutoff frequency
% b = fir1(50,f1/fN); % lowpass Hamming window filter (default)
% ufilt = nan(size(u));
% vfilt = nan(size(u));
% for ii = 1:length(depth_MOOR)
%     ufilt(ii,:) = filtfilt(b,1,u(ii,:));
%     vfilt(ii,:) = filtfilt(b,1,v(ii,:));
% end
r = round(.5/dt);
dt = r*dt;
for ii = 1:length(depth_MOOR)
    ufilt(ii,:) = decimate(u(ii,:),r);
    vfilt(ii,:) = decimate(v(ii,:),r);
end

u = ufilt;
v = vfilt;

%% load mooring data: WOCE ACM1 ACCP2 312
% (02 Oct 93 - 29 Oct 95)

% WOCE ACM1 mooring data is here
MOORdir = '~worthamc/Documents/research/data/currentmetersOSU/Atlantic/WOCE/acm1/accp2/';

filenames = {'rcm02736.str' 'rcm02738.str' 'rcm02742.str' 'rcm02744.str' 'rcm02746.str' 'rcm02748.str'};

% get the data
for ii = 1:length(filenames);
    fid = fopen([MOORdir filenames{ii}],'r');
    % collect information from header
    headlines = textscan(fid,'%d header lines',1);headlines = headlines{:};
    datalines = textscan(fid,'%d data lines',1);datalines = datalines{:};
    textscan(fid,'%s',1); % skipping a line
    expname = textscan(fid,'Experiment name: %[^\n]',1);expname = char(deblank(expname{1}));
    moorname = textscan(fid,'Mooring name: %[^\n]',1);moorname = char(deblank(moorname{1}));
    latlong = textscan(fid,'Mooring position: %f deg %c, %f deg %c',1);
    if (latlong{2}=='N'), lat = latlong{1}; else lat = -latlong{1}; end
    if (latlong{4}=='E'), lon = latlong{3}; else lon = 360-latlong{3}; end
    instdepth = textscan(fid,'Instrument depth: %d meters',1);instdepth = instdepth{1};
    waterdepth = textscan(fid,'Seafloor depth: %d meters',1);waterdepth = waterdepth{1};
    instrument = textscan(fid,'Instrument type: %[^\n]',1);instrument = char(deblank(instrument{1}));
    accession = textscan(fid,'CMDB accession number: %d',1);accession = accession{1};
    for i=1:5
        textscan(fid,'%s',1); % skipping some lines
    end
    velunits = textscan(fid,'speed (%[^)]',1);velunits = char(velunits{1});
    % load the data
    fseek(fid,0,'bof'); %reset file position indicator to beginning
    if headlines==20 % setting number of fields based on header length
        fields='%f %f %f %f %f %f %f %f %f';
    elseif headlines==21
        fields='%f %f %f %f %f %f %f %f %f %f';
    elseif headlines==22
        fields='%f %f %f %f %f %f %f %f %f %f %f';
    elseif headlines==23
        fields='%f %f %f %f %f %f %f %f %f %f %f %f';
    else
        disp('check number of fields')
    end
    theData = textscan(fid, fields, datalines, 'headerlines', headlines);
    fclose(fid);
    
    % add to final data vectors
    depth_MOOR(ii) = -double(instdepth);
    if datalines==1515;
        startind = 2;
    else
        startind = 1;
    end    
    u(ii,:) = theData{7}(startind:end);
    v(ii,:) = theData{8}(startind:end);
    hr = theData{1}(startind:end);
    mm = mod(hr,100);
    hr = (hr-mm)/100;
    time = datenum(theData{4}(startind:end)+1900,theData{3}(startind:end),theData{2}(startind:end),hr,mm,0);
end

dt = (time(2) - time(1)); % days

% remove time mean velocity
ubar = mean(u,2); ubar = repmat(ubar,[1,length(time)]);
vbar = mean(v,2); vbar = repmat(vbar,[1,length(time)]);
u = u - ubar;
v = v - vbar;

% make units m/s
u = u/100;
v = v/100;

% fudge latitude just a little north...mooring is too close to Abacco
% Island for interpolation from OCCA dataset to work.
lat = lat+.1;

% no filtering... time step is already 12 hours

%% load mooring data: WOCE ACM1 WATTS M271
% (22 Jun 90 - 09 Feb 92)
% removed instrument 6... stops after 6 months

% WOCE ACM1 mooring data is here
MOORdir = '~worthamc/Documents/research/data/currentmetersOSU/Atlantic/WOCE/acm1/watts/';

filenames = {'rcm02531.str' 'rcm02533.str' 'rcm02534.str' 'rcm02535.str' 'rcm02537.str' 'rcm02539.str'};

% get the data
for ii = 1:length(filenames);
    fid = fopen([MOORdir filenames{ii}],'r');
    % collect information from header
    headlines = textscan(fid,'%d header lines',1);headlines = headlines{:};
    datalines = textscan(fid,'%d data lines',1);datalines = datalines{:};
    textscan(fid,'%s',1); % skipping a line
    expname = textscan(fid,'Experiment name: %[^\n]',1);expname = char(deblank(expname{1}));
    moorname = textscan(fid,'Mooring name: %[^\n]',1);moorname = char(deblank(moorname{1}));
    latlong = textscan(fid,'Mooring position: %f deg %c, %f deg %c',1);
    if (latlong{2}=='N'), lat = latlong{1}; else lat = -latlong{1}; end
    if (latlong{4}=='E'), lon = latlong{3}; else lon = 360-latlong{3}; end
    instdepth = textscan(fid,'Instrument depth: %d meters',1);instdepth = instdepth{1};
    waterdepth = textscan(fid,'Seafloor depth: %d meters',1);waterdepth = waterdepth{1};
    instrument = textscan(fid,'Instrument type: %[^\n]',1);instrument = char(deblank(instrument{1}));
    accession = textscan(fid,'CMDB accession number: %d',1);accession = accession{1};
    for i=1:5
        textscan(fid,'%s',1); % skipping some lines
    end
    velunits = textscan(fid,'speed (%[^)]',1);velunits = char(velunits{1});
    % load the data
    fseek(fid,0,'bof'); %reset file position indicator to beginning
    if headlines==20 % setting number of fields based on header length
        fields='%f %f %f %f %f %f %f %f %f';
    elseif headlines==21
        fields='%f %f %f %f %f %f %f %f %f %f';
    elseif headlines==22
        fields='%f %f %f %f %f %f %f %f %f %f %f';
    elseif headlines==23
        fields='%f %f %f %f %f %f %f %f %f %f %f %f';
    else
        disp('check number of fields')
    end
    theData = textscan(fid, fields, datalines, 'headerlines', headlines);
    fclose(fid);
    
    % add to final data vectors
    depth_MOOR(ii) = -double(instdepth);
    stopind = 1194;
    u(ii,:) = theData{7}(1:stopind);
    v(ii,:) = theData{8}(1:stopind);
    hr = theData{1}(1:stopind);
    mm = mod(hr,100);
    hr = (hr-mm)/100;
    time = datenum(theData{4}(1:stopind)+1900,theData{3}(1:stopind),theData{2}(1:stopind),hr,mm,0);
end

dt = (time(2) - time(1)); % days

% remove time mean velocity
ubar = mean(u,2); ubar = repmat(ubar,[1,length(time)]);
vbar = mean(v,2); vbar = repmat(vbar,[1,length(time)]);
u = u - ubar;
v = v - vbar;

% make units m/s
u = u/100;
v = v/100;

% no filtering... time step is already 12 hours

%% load mooring data: LLWODP-W CMMW-10
% (22 Sep 82 - 02 Sep 83)

% LLWODP-w mooring data is here
MOORdir = '~worthamc/Documents/research/data/currentmetersOSU/Pacific/llwodp_w/';

filenames = {'rcm00198.str' 'rcm00199.str' 'rcm00200.str' 'rcm00201.str' 'rcm00202.str' 'rcm00203.str'};

% get the data
for ii = 1:length(filenames);
    fid = fopen([MOORdir filenames{ii}],'r');
    % collect information from header
    headlines = textscan(fid,'%d header lines',1);headlines = headlines{:};
    datalines = textscan(fid,'%d data lines',1);datalines = datalines{:};
    textscan(fid,'%s',1); % skipping a line
    expname = textscan(fid,'Experiment name: %[^\n]',1);expname = char(deblank(expname{1}));
    moorname = textscan(fid,'Mooring name: %[^\n]',1);moorname = char(deblank(moorname{1}));
    latlong = textscan(fid,'Mooring position: %f deg %c, %f deg %c',1);
    if (latlong{2}=='N'), lat = latlong{1}; else lat = -latlong{1}; end
    if (latlong{4}=='E'), lon = latlong{3}; else lon = 360-latlong{3}; end
    instdepth = textscan(fid,'Instrument depth: %d meters',1);instdepth = instdepth{1};
    waterdepth = textscan(fid,'Seafloor depth: %d meters',1);waterdepth = waterdepth{1};
    instrument = textscan(fid,'Instrument type: %[^\n]',1);instrument = char(deblank(instrument{1}));
    accession = textscan(fid,'CMDB accession number: %d',1);accession = accession{1};
    for i=1:5
        textscan(fid,'%s',1); % skipping some lines
    end
    velunits = textscan(fid,'speed (%[^)]',1);velunits = char(velunits{1});
    % load the data
    fseek(fid,0,'bof'); %reset file position indicator to beginning
    if headlines==20 % setting number of fields based on header length
        fields='%f %f %f %f %f %f %f %f %f';
    elseif headlines==21
        fields='%f %f %f %f %f %f %f %f %f %f';
    elseif headlines==22
        fields='%f %f %f %f %f %f %f %f %f %f %f';
    elseif headlines==23
        fields='%f %f %f %f %f %f %f %f %f %f %f %f';
    else
        disp('check number of fields')
    end
    theData = textscan(fid, fields, datalines, 'headerlines', headlines);
    fclose(fid);
    
    % interpolate velocity for the mooring taking data every 30min
    if ~(datalines==8297)
        hr = theData{1};
        mm = mod(hr,100);
        hr = (hr-mm)/100;
        time_30min = datenum(theData{4}+1900,theData{3},theData{2},hr,mm,0);
        theData{7} = interp1(time_30min,theData{7},time,'linear','extrap');
        theData{8} = interp1(time_30min,theData{8},time,'linear','extrap');
    end
    
    % add to final data vectors
    depth_MOOR(ii) = -double(instdepth);
    stopind = 1194;
    u(ii,:) = theData{7};
    v(ii,:) = theData{8};
    hr = theData{1};
    mm = mod(hr,100);
    hr = (hr-mm)/100;
    time = datenum(theData{4}+1900,theData{3},theData{2},hr,mm,0);
end

dt = (time(2) - time(1)); % days

% remove time mean velocity
ubar = mean(u,2); ubar = repmat(ubar,[1,length(time)]);
vbar = mean(v,2); vbar = repmat(vbar,[1,length(time)]);
u = u - ubar;
v = v - vbar;

% make units m/s
u = u/100;
v = v/100;

% lowpass filter
% fN = 1/2/dt; % Nyquist frequency (cph)
% f1 = 1/24; % cutoff frequency
% b = fir1(50,f1/fN); % lowpass Hamming window filter (default)
% ufilt = nan(size(u));
% vfilt = nan(size(u));
% for ii = 1:length(depth_MOOR)
%     ufilt(ii,:) = filtfilt(b,1,u(ii,:));
%     vfilt(ii,:) = filtfilt(b,1,v(ii,:));
% end
% filter and decimate to twice daily
r = round(.5/dt);
dt = r*dt;
for ii = 1:length(depth_MOOR)
    ufilt(ii,:) = decimate(u(ii,:),r);
    vfilt(ii,:) = decimate(v(ii,:),r);
end

u = ufilt;
v = vfilt;

%% load mooring data: pacmisc WHOI 793
% (26 Oct 83 - 30 Sep 85)

% pacmisc mooring data is here
MOORdir = '~worthamc/Documents/research/data/currentmetersOSU/Pacific/pacmisc/';

filenames = {'rcm03221.str' 'rcm03222.str' 'rcm03223.str'};

% get the data
for ii = 1:length(filenames);
    fid = fopen([MOORdir filenames{ii}],'r');
    % collect information from header
    headlines = textscan(fid,'%d header lines',1);headlines = headlines{:};
    datalines = textscan(fid,'%d data lines',1);datalines = datalines{:};
    textscan(fid,'%s',1); % skipping a line
    expname = textscan(fid,'Experiment name: %[^\n]',1);expname = char(deblank(expname{1}));
    moorname = textscan(fid,'Mooring name: %[^\n]',1);moorname = char(deblank(moorname{1}));
    latlong = textscan(fid,'Mooring position: %f deg %c, %f deg %c',1);
    if (latlong{2}=='N'), lat = latlong{1}; else lat = -latlong{1}; end
    if (latlong{4}=='E'), lon = latlong{3}; else lon = 360-latlong{3}; end
    instdepth = textscan(fid,'Instrument depth: %d meters',1);instdepth = instdepth{1};
    waterdepth = textscan(fid,'Seafloor depth: %d meters',1);waterdepth = waterdepth{1};
    instrument = textscan(fid,'Instrument type: %[^\n]',1);instrument = char(deblank(instrument{1}));
    accession = textscan(fid,'CMDB accession number: %d',1);accession = accession{1};
    for i=1:5
        textscan(fid,'%s',1); % skipping some lines
    end
    velunits = textscan(fid,'speed (%[^)]',1);velunits = char(velunits{1});
    % load the data
    fseek(fid,0,'bof'); %reset file position indicator to beginning
    if headlines==20 % setting number of fields based on header length
        fields='%f %f %f %f %f %f %f %f %f';
    elseif headlines==21
        fields='%f %f %f %f %f %f %f %f %f %f';
    elseif headlines==22
        fields='%f %f %f %f %f %f %f %f %f %f %f';
    elseif headlines==23
        fields='%f %f %f %f %f %f %f %f %f %f %f %f';
    else
        disp('check number of fields')
    end
    theData = textscan(fid, fields, datalines, 'headerlines', headlines);
    fclose(fid);
    
    if ~(datalines==33841)
        theData{1} = theData{1}(1:33841);
        theData{2} = theData{2}(1:33841);
        theData{3} = theData{3}(1:33841);
        theData{4} = theData{4}(1:33841);
        theData{7} = theData{7}(1:33841);
        theData{8} = theData{8}(1:33841);
    end
    
    % add to final data vectors
    depth_MOOR(ii) = -double(instdepth);
    u(ii,:) = theData{7};
    v(ii,:) = theData{8};
    hr = theData{1};
    mm = mod(hr,100);
    hr = (hr-mm)/100;
    time = datenum(theData{4}+1900,theData{3},theData{2},hr,mm,0);
end

dt = (time(2) - time(1)); % days

% remove time mean velocity
ubar = mean(u,2); ubar = repmat(ubar,[1,length(time)]);
vbar = mean(v,2); vbar = repmat(vbar,[1,length(time)]);
u = u - ubar;
v = v - vbar;

% make units m/s
u = u/100;
v = v/100;

% lowpass filter
% fN = 1/2/dt; % Nyquist frequency (cph)
% f1 = 1/24; % cutoff frequency
% b = fir1(50,f1/fN); % lowpass Hamming window filter (default)
% ufilt = nan(size(u));
% vfilt = nan(size(u));
% for ii = 1:length(depth_MOOR)
%     ufilt(ii,:) = filtfilt(b,1,u(ii,:));
%     vfilt(ii,:) = filtfilt(b,1,v(ii,:));
% end
% filter and decimate to twice daily
r = round(.5/dt);
dt = r*dt;
for ii = 1:length(depth_MOOR)
    ufilt(ii,:) = decimate(u(ii,:),r);
    vfilt(ii,:) = decimate(v(ii,:),r);
end

u = ufilt;
v = vfilt;

%% load mooring data: vents 88V22
% (14 Jul 88 - 13 Aug 89)

% vents mooring data is here
MOORdir = '~worthamc/Documents/research/data/currentmetersOSU/Pacific/vents/';

filenames = {'acc02984.str' 'acc02985.str' 'acc02986.str' 'acc02987.str' 'acc02988.str' 'acc02989.str'};

% get the data
for ii = 1:length(filenames);
    fid = fopen([MOORdir filenames{ii}],'r');
    % collect information from header
    headlines = textscan(fid,'%d header lines',1);headlines = headlines{:};
    datalines = textscan(fid,'%d data lines',1);datalines = datalines{:};
    textscan(fid,'%s',1); % skipping a line
    expname = textscan(fid,'Experiment name: %[^\n]',1);expname = char(deblank(expname{1}));
    moorname = textscan(fid,'Mooring name: %[^\n]',1);moorname = char(deblank(moorname{1}));
    latlong = textscan(fid,'Mooring position: %f deg %c, %f deg %c',1);
    if (latlong{2}=='N'), lat = latlong{1}; else lat = -latlong{1}; end
    if (latlong{4}=='E'), lon = latlong{3}; else lon = 360-latlong{3}; end
    instdepth = textscan(fid,'Instrument depth: %d meters',1);instdepth = instdepth{1};
    waterdepth = textscan(fid,'Seafloor depth: %d meters',1);waterdepth = waterdepth{1};
    instrument = textscan(fid,'Instrument type: %[^\n]',1);instrument = char(deblank(instrument{1}));
    accession = textscan(fid,'CMDB accession number: %d',1);accession = accession{1};
    for i=1:5
        textscan(fid,'%s',1); % skipping some lines
    end
    velunits = textscan(fid,'speed (%[^)]',1);velunits = char(velunits{1});
    % load the data
    fseek(fid,0,'bof'); %reset file position indicator to beginning
    if headlines==20 % setting number of fields based on header length
        fields='%f %f %f %f %f %f %f %f %f';
    elseif headlines==21
        fields='%f %f %f %f %f %f %f %f %f %f';
    elseif headlines==22
        fields='%f %f %f %f %f %f %f %f %f %f %f';
    elseif headlines==23
        fields='%f %f %f %f %f %f %f %f %f %f %f %f';
    else
        disp('check number of fields')
    end
    theData = textscan(fid, fields, datalines, 'headerlines', headlines);
    fclose(fid);
    
    % add to final data vectors
    depth_MOOR(ii) = -double(instdepth);
    u(ii,:) = theData{7};
    v(ii,:) = theData{8};
    hr = theData{1};
    mm = mod(hr,100);
    hr = (hr-mm)/100;
    time = datenum(theData{4}+1900,theData{3},theData{2},hr,mm,0);
end

dt = (time(2) - time(1)); % days

% remove time mean velocity
ubar = mean(u,2); ubar = repmat(ubar,[1,length(time)]);
vbar = mean(v,2); vbar = repmat(vbar,[1,length(time)]);
u = u - ubar;
v = v - vbar;

% make units m/s
u = u/100;
v = v/100;

% lowpass filter
% fN = 1/2/dt; % Nyquist frequency (cph)
% f1 = 1/24; % cutoff frequency
% b = fir1(50,f1/fN); % lowpass Hamming window filter (default)
% ufilt = nan(size(u));
% vfilt = nan(size(u));
% for ii = 1:length(depth_MOOR)
%     ufilt(ii,:) = filtfilt(b,1,u(ii,:));
%     vfilt(ii,:) = filtfilt(b,1,v(ii,:));
% end
% filter and decimate to twice daily
r = round(.5/dt);
dt = r*dt;
for ii = 1:length(depth_MOOR)
    ufilt(ii,:) = decimate(u(ii,:),r);
    vfilt(ii,:) = decimate(v(ii,:),r);
end

u = ufilt;
v = vfilt;

%% load mooring data: vents 94V86
% (16 Sep 94 - 17 Jun 95)

% vents mooring data is here
MOORdir = '~worthamc/Documents/research/data/currentmetersOSU/Pacific/vents/';

filenames = {'acc05417.str' 'acc05418.str' 'acc05419.str' 'acc05420.str'};

% get the data
for ii = 1:length(filenames);
    fid = fopen([MOORdir filenames{ii}],'r');
    % collect information from header
    headlines = textscan(fid,'%d header lines',1);headlines = headlines{:};
    datalines = textscan(fid,'%d data lines',1);datalines = datalines{:};
    textscan(fid,'%s',1); % skipping a line
    expname = textscan(fid,'Experiment name: %[^\n]',1);expname = char(deblank(expname{1}));
    moorname = textscan(fid,'Mooring name: %[^\n]',1);moorname = char(deblank(moorname{1}));
    latlong = textscan(fid,'Mooring position: %f deg %c, %f deg %c',1);
    if (latlong{2}=='N'), lat = latlong{1}; else lat = -latlong{1}; end
    if (latlong{4}=='E'), lon = latlong{3}; else lon = 360-latlong{3}; end
    instdepth = textscan(fid,'Instrument depth: %d meters',1);instdepth = instdepth{1};
    waterdepth = textscan(fid,'Seafloor depth: %d meters',1);waterdepth = waterdepth{1};
    instrument = textscan(fid,'Instrument type: %[^\n]',1);instrument = char(deblank(instrument{1}));
    accession = textscan(fid,'CMDB accession number: %d',1);accession = accession{1};
    for i=1:5
        textscan(fid,'%s',1); % skipping some lines
    end
    velunits = textscan(fid,'speed (%[^)]',1);velunits = char(velunits{1});
    % load the data
    fseek(fid,0,'bof'); %reset file position indicator to beginning
    if headlines==20 % setting number of fields based on header length
        fields='%f %f %f %f %f %f %f %f %f';
    elseif headlines==21
        fields='%f %f %f %f %f %f %f %f %f %f';
    elseif headlines==22
        fields='%f %f %f %f %f %f %f %f %f %f %f';
    elseif headlines==23
        fields='%f %f %f %f %f %f %f %f %f %f %f %f';
    else
        disp('check number of fields')
    end
    theData = textscan(fid, fields, datalines, 'headerlines', headlines);
    fclose(fid);
    
    % add to final data vectors
    depth_MOOR(ii) = -double(instdepth);
    u(ii,:) = theData{7};
    v(ii,:) = theData{8};
    hr = theData{1};
    mm = mod(hr,100);
    hr = (hr-mm)/100;
    time = datenum(theData{4}+1900,theData{3},theData{2},hr,mm,0);
end

dt = (time(2) - time(1)); % days

% remove time mean velocity
ubar = mean(u,2); ubar = repmat(ubar,[1,length(time)]);
vbar = mean(v,2); vbar = repmat(vbar,[1,length(time)]);
u = u - ubar;
v = v - vbar;

% make units m/s
u = u/100;
v = v/100;

% lowpass filter
% fN = 1/2/dt; % Nyquist frequency (cph)
% f1 = 1/24; % cutoff frequency
% b = fir1(50,f1/fN); % lowpass Hamming window filter (default)
% ufilt = nan(size(u));
% vfilt = nan(size(u));
% for ii = 1:length(depth_MOOR)
%     ufilt(ii,:) = filtfilt(b,1,u(ii,:));
%     vfilt(ii,:) = filtfilt(b,1,v(ii,:));
% end
% filter and decimate to twice daily
r = round(.5/dt);
dt = r*dt;
for ii = 1:length(depth_MOOR)
    ufilt(ii,:) = decimate(u(ii,:),r);
    vfilt(ii,:) = decimate(v(ii,:),r);
end

u = ufilt;
v = vfilt;

%% load mooring data: WestPac WHOI 704
% (19 Jul 80 - 12 May 81)

% westpac mooring data is here
MOORdir = '~worthamc/Documents/research/data/currentmetersOSU/Pacific/wp/';

filenames = {'rcm01292.str' 'rcm01293.str' 'rcm01294.str' 'rcm01295.str' 'rcm01296.str'};

% get the data
for ii = 1:length(filenames);
    fid = fopen([MOORdir filenames{ii}],'r');
    % collect information from header
    headlines = textscan(fid,'%d header lines',1);headlines = headlines{:};
    datalines = textscan(fid,'%d data lines',1);datalines = datalines{:};
    textscan(fid,'%s',1); % skipping a line
    expname = textscan(fid,'Experiment name: %[^\n]',1);expname = char(deblank(expname{1}));
    moorname = textscan(fid,'Mooring name: %[^\n]',1);moorname = char(deblank(moorname{1}));
    latlong = textscan(fid,'Mooring position: %f deg %c, %f deg %c',1);
    if (latlong{2}=='N'), lat = latlong{1}; else lat = -latlong{1}; end
    if (latlong{4}=='E'), lon = latlong{3}; else lon = 360-latlong{3}; end
    instdepth = textscan(fid,'Instrument depth: %d meters',1);instdepth = instdepth{1};
    waterdepth = textscan(fid,'Seafloor depth: %d meters',1);waterdepth = waterdepth{1};
    instrument = textscan(fid,'Instrument type: %[^\n]',1);instrument = char(deblank(instrument{1}));
    accession = textscan(fid,'CMDB accession number: %d',1);accession = accession{1};
    for i=1:5
        textscan(fid,'%s',1); % skipping some lines
    end
    velunits = textscan(fid,'speed (%[^)]',1);velunits = char(velunits{1});
    % load the data
    fseek(fid,0,'bof'); %reset file position indicator to beginning
    if headlines==20 % setting number of fields based on header length
        fields='%f %f %f %f %f %f %f %f %f';
    elseif headlines==21
        fields='%f %f %f %f %f %f %f %f %f %f';
    elseif headlines==22
        fields='%f %f %f %f %f %f %f %f %f %f %f';
    elseif headlines==23
        fields='%f %f %f %f %f %f %f %f %f %f %f %f';
    else
        disp('check number of fields')
    end
    theData = textscan(fid, fields, datalines, 'headerlines', headlines);
    fclose(fid);
    
    % add to final data vectors
    depth_MOOR(ii) = -double(instdepth);
    u(ii,:) = theData{7};
    v(ii,:) = theData{8};
    hr = theData{1};
    mm = mod(hr,100);
    hr = (hr-mm)/100;
    time = datenum(theData{4}+1900,theData{3},theData{2},hr,mm,0);
end

dt = (time(2) - time(1)); % days

% remove time mean velocity
ubar = mean(u,2); ubar = repmat(ubar,[1,length(time)]);
vbar = mean(v,2); vbar = repmat(vbar,[1,length(time)]);
u = u - ubar;
v = v - vbar;

% make units m/s
u = u/100;
v = v/100;

% lowpass filter
% fN = 1/2/dt; % Nyquist frequency (cph)
% f1 = 1/24; % cutoff frequency
% b = fir1(50,f1/fN); % lowpass Hamming window filter (default)
% ufilt = nan(size(u));
% vfilt = nan(size(u));
% for ii = 1:length(depth_MOOR)
%     ufilt(ii,:) = filtfilt(b,1,u(ii,:));
%     vfilt(ii,:) = filtfilt(b,1,v(ii,:));
% end
% filter and decimate to twice daily
r = round(.5/dt);
dt = r*dt;
for ii = 1:length(depth_MOOR)
    ufilt(ii,:) = decimate(u(ii,:),r);
    vfilt(ii,:) = decimate(v(ii,:),r);
end

u = ufilt;
v = vfilt;

%% get OCCA data

% load OCCA grid (check locations!)
lat_OCCA = nc_varget([OCCAdir 'DDtheta.0406annclim.nc'],'Latitude_t');
lon_OCCA = nc_varget([OCCAdir 'DDtheta.0406annclim.nc'],'Longitude_t');
depth_OCCA = -nc_varget([OCCAdir 'DDtheta.0406annclim.nc'],'Depth_c')';

% load OCCA data
salt_OCCA = nc_varget([OCCAdir 'DDsalt.0406annclim.nc'],'salt');
theta_OCCA = nc_varget([OCCAdir 'DDtheta.0406annclim.nc'],'theta');

% interpolate OCCA data to mooring location
S = zeros(1,length(depth_OCCA));
theta = zeros(1,length(depth_OCCA));
for ii=1:length(depth_OCCA)
    S(ii) = interp2(lon_OCCA,lat_OCCA,squeeze(salt_OCCA(ii,:,:)),lon,lat);
    theta(ii) = interp2(lon_OCCA,lat_OCCA,squeeze(theta_OCCA(ii,:,:)),lon,lat);
end

% clear OCCA data
clear salt_OCCA theta_OCCA lat_OCCA lon_OCCA

% cut off nans
depth_OCCA = depth_OCCA(~isnan(theta));
S = S(~isnan(theta));
theta = theta(~isnan(theta));

%% get some other climatology... does it make a difference?

%% compute BTT mode shapes

% convert to in situ temp
pres = sw_pres(-depth_OCCA,lat);
T = sw_temp(S,theta,pres,0); clear pres;

% compute BTT mode shapes
iplot = 1;
nplot = 5;
[Ld,eigenvecs,rho,N2cph] = BTTmodes(S,T,depth_OCCA,lat,iplot,nplot);

% interpolate mode shapes to current meter depths
eigenvecs_int = interp1(depth_OCCA,eigenvecs,depth_MOOR,'linear','extrap');
nmodes = 5;

% plot
figure(99);
subplot(1,4,1:2)
hold on
plot(eigenvecs_int(:,1:nmodes),depth_MOOR,'--')
plot(0,depth_MOOR,'ko')
hold off
title(moorname)

%% modal decomposition 1: Least Squares

% E_temp = eigenvecs_int(:,1:nmodes);
% alpha_u = nan(nmodes,length(time));
% alpha_v = nan(nmodes,length(time));
% for tt = 1:length(time)
%     
%     % u velocity
%     u_temp = u(:,tt);
%     alpha_u(:,tt) = E_temp\u_temp;
%     
%     % v velocity
%     v_temp = v(:,tt);
%     alpha_v(:,tt) = E_temp\v_temp;
% end
% vectorize!  mldivide does it right
alphaLS_u = eigenvecs_int(:,1:nmodes)\u;
alphaLS_v = eigenvecs_int(:,1:nmodes)\v;

%% modal decomposition 2: projection

dz = diff(depth_MOOR);
dz = [depth_MOOR(1) dz];
dz = repmat(dz',[1 5]);
alphaP_u = 1/depth_MOOR(end)*(eigenvecs_int(:,1:5).*dz)'*u;
alphaP_v = 1/depth_MOOR(end)*(eigenvecs_int(:,1:5).*dz)'*v;

% dz = diff(depth_MOOR);
% dz = [depth_MOOR(1) dz];
% alphaP_u = nan(nmodes,length(time));
% alphaP_v = nan(nmodes,length(time));
% for nn = 1:nmodes
%     E_temp = repmat(eigenvecs_int(:,nn),[1 length(time)]);
%     alphaP_u(nn,:) = 1/depth_MOOR(end)*trapz(dz,u.*E_temp);
%     alphaP_v(nn,:) = 1/depth_MOOR(end)*trapz(dz,v.*E_temp);
% end

%% modal decomposition 3: Gauss-Markov estimate

E2u = mean(var(u,0,2));
E2v = mean(var(v,0,2));
P0u = E2u * diag([1 1 1/2 1/4 1/8]) /2.88;
P0v = E2v * diag([1 1 1/2 1/4 1/8]) /2.88;
A = eigenvecs_int(:,1:5);
sigma2 = 9e-4;

alphaGM_u = P0u*A' * inv(A*P0u*A' + sigma2*eye(size(A*P0u*A'))) * u;
alphaGM_v = P0v*A' * inv(A*P0v*A' + sigma2*eye(size(A*P0v*A'))) * v;

%% analysis

% choose method to analyze
% alpha_u = alphaLS_u;
% alpha_v = alphaLS_v;

alpha_u = alphaGM_u;
alpha_v = alphaGM_v;

% root mean (space and time) square for each mode, normalized
alpha_uRMS = zeros(nmodes,1);
alpha_vRMS = zeros(nmodes,1);
for ii = 1:nmodes
    alpha_uRMS(ii) = sqrt( nanmean( alpha_u(ii,:).^2 ) );
    alpha_vRMS(ii) = sqrt( nanmean( alpha_v(ii,:).^2 ) );
end
alpha_uRMS = alpha_uRMS./max(alpha_uRMS);
alpha_uRMS
alpha_vRMS = alpha_vRMS./max(alpha_vRMS);
alpha_vRMS

% reconstructed velocities
tt = 100;
fno = fno+1;
figure(fno);clf
subplot(121)
plot(u(:,tt),depth_MOOR,v(:,tt),depth_MOOR)
title(['data, ' moorname])
xlabel('velocity (m/s)')
ylabel('depth (m)')
subplot(122)
plot(eigenvecs_int(:,1:nmodes)*alpha_u(:,tt),depth_MOOR,eigenvecs_int(:,1:nmodes)*alpha_v(:,tt),depth_MOOR)
hold on
plot(u(:,tt)-eigenvecs_int(:,1:nmodes)*alpha_u(:,tt),depth_MOOR,'--',v(:,tt)-eigenvecs_int(:,1:nmodes)*alpha_v(:,tt),depth_MOOR,'--')
hold off
title('reconstructed')
xlabel('velocity (m/s)')

% variance captured by first few modes
a0VAR = var( alpha_u(1,:)'*eigenvecs_int(:,1)' ); % variance in BT mode, assuming no coupling
a1VAR = var( alpha_u(2,:)'*eigenvecs_int(:,2)' ); % variance in BC1 mode, assuming no coupling
aVAR = var( alpha_u'*eigenvecs_int(:,1:nmodes)' ); % variance in BT-BC4 modes combined
uVAR = var( u,0,2); % variance in zonal velocity at this location
fno = fno+1;
figure(fno);clf;
plot(uVAR,depth_MOOR,'k')
hold on
plot(aVAR,depth_MOOR,'b')
%plot(a0VAR,depth_MOOR,'b:')
%plot(a1VAR,depth_MOOR,'b--')
hold off
legend('current meters','BT-BC4','BT','BC1','location','southeast')
xlabel('variance (m^2/s^2)')
ylabel('depth (m)')
% annotate
text(max(uVAR)*.7,min(depth_MOOR),num2str(alpha_uRMS,2),'HorizontalAlignment','left','fontsize',12)
text(max(uVAR)*.7,min(depth_MOOR)*.8,'\alpha_{u,RMS}','horizontalalignment','left','fontsize',12)
title(moorname)

% modal coupling
[r,p] = corrcoef( alpha_u(1,:),alpha_u(2,:));

% modal coherence
fno = fno+1;
figure(fno);clf
NW=8;qbias=1;confn=0;qplot=1;
[s,c,ph,ci,phi] = cmtm( alpha_u(1,:),alpha_u(2,:),dt,NW,qbias,confn,qplot);

%% map

mooringnames = {'lotus 1' 'lotus 2' 'ACM1 ACCP2 312' 'ACM1 WATTS M271' 'LLWODP-W CMMW-10' 'WHOI 793' 'vents 88V22' 'vents 94V86' 'WHOI 704'}';
latitudes = [33.9000 33.9000 26.4980 26.5000 39.4730 38.9600 45.4560 44.9800 27.9930]';
longitudes = [290.0000 290.0000 284.3170 284.3170 232.3180 207.9700 228.8230 229.8270 151.9400]';

fno = fno+1;
figure(fno);clf
m_proj('Miller Cylindrical','lon',[0 360],'lat',[-66 66]);
m_plot(longitudes,latitudes,'ko','markerfacecolor','k')
set(findobj('tag','m_grid_color'),'facecolor','none'); % necessary for m_countourf plots
title('map of moorings')
% add bathymetry
[cs,h] = m_etopo2('contour',[-4000 -2000],'edgecolor',[.7 .7 .7]);
m_coast('patch',[.5 .5 .5],'edgecolor',[.5 .5 .5]);
m_grid('box','on','tickdir','in');
m_text(longitudes+1.5,latitudes-2,mooringnames)





