function E = convert_to_SRF(data, index)
%CONVERT_TO_SRF Summary of this function goes here
%   Detailed explanation goes here
if ~exist('index', 'var')
    index = 1;
end
if ~exist('tofile', 'var')
    tofile = 0;
end
ep = data.epoch;
time = ep(index);
tds_mode = convertCharsToStrings(data.tds_config_label(:,index));
% Pachenko's antenna angle
pacang = 132.2;
% RPW ANT in SRF 
V1=[0, 1];
% Geometric antenna
%V2=[0.81915204586197221, -0.57357643410459713]; 
%V3=[-0.81915204586197221, -0.57357643410459713]; 
% Pachenko antenna
V2=[sind(pacang), cosd(pacang)]; 
V3=[-sind(pacang), cosd(pacang)]; 
V2 = V2./norm(V2);
V3 = V3./norm(V3);
% XLD1 TDS mode
V1V3 = V1 - V3; V1V3 = V1V3/ norm(V1V3); %CH1
V2V1 = V2 - V1; V2V1 = V2V1/ norm(V2V1); %CH2
if contains(tds_mode, 'SE1')
% SE1 TDS mode
    M(1,:) = V1; %CH1
    M(2,:) = V2; %CH2
else 
    % construct transformation matrix ANT->SRF 
    M(1,:) = V1V3; % V1-V3
    M(2,:) = V2V1; % V2-V1
    % DIFF1 TDS mode is same
end
time = ep(index);
% load TDS TSWF data
% if it's at the OKF then use find_file
nsamp = double(data.samples_per_ch(index));
fs = data.samp_rate(index);
tvec = (0:nsamp-1) / fs * 1e3;
ww = data.data(:,1:nsamp,index);
% projection: E = MAT(ANT->SRF) * V; where MAT(2,2) and V is observed field 
M = inv(M);
E = M*ww(1:2,:) * 1e3; % transformation into SRF (Y-Z)
end

