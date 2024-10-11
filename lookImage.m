clear
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%% Params settings START %%%%%%%%%%%%%%%%%%%%%%%%%%%
myname = 'infi1';
filename = sprintf('%s_H.mci',myname);
disp(['loading ' filename])
fid = fopen(filename, 'r');
A = fscanf(fid,'%f',[1 Inf])';
fclose(fid);
n = 1;
Nphotons = A(n); n = n + 1;
p = A(n); n = n + 1;
Ndetectors = A(n); n = n + 1;
det_radius = A(n); n = n + 1;
cos_accept = A(n); n = n + 1;
Nx = A(n); n = n + 1;
Ny = A(n); n = n + 1;
Nz = A(n); n = n + 1;
dx = A(n); n = n + 1;
dy = A(n); n = n + 1;
dz = A(n); n = n + 1;
radius = A(n); n = n + 1;
zsurf = A(n); n = n + 1;
Nt = A(n);  n = n + 1;
lambda_start = 1250e-7;% cm, start wavelength, make sure it is the same as GenTissue.m
lambda_stop = 1350e-7;% cm, end wavelength, make sure it is the same as GenTissue.m
samplePoints = 1024;% sampling of wavelength, make sure it is the same as GenTissue.m
nm = 1.33; % Refractive index, make sure it is the same as GenTissue.m
lambda_lin = linspace(lambda_start,lambda_stop,samplePoints);
k_nolin = 2*pi./lambda_lin;
k_lin = linspace(2*pi/lambda_start,2*pi/lambda_stop,samplePoints);
Bscan_number = 4; % make sure it is the same as fullwaveMC_Bscan_OCTA.c
sig_type = 'sig'; %sig or sig_ss or sig_ms
noise_amp = 0.05; % noise amplitude
kernel_size = 7; % OCTA window size
u1 = 0.5;%mm/s, averaged velocity, make sure it is the same as Gen_particles.m
numBins = 20; % Q-OCTA binning number
lower_bound_num = 25; % #%,Q-OCTA variance range
upper_bound_num = 75; % #%,Q-OCTA variance range
velocity_threshold = 4.5;% mm/s, threshold to reduce the variance
decorrelation_threshold = 0.7;% 1, threshold to reduce the variance
load('parameters/vessel_information.mat')
%% %%%%%%%%%%%%%%%%%%%%% Params settings DONE %%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OCT Bscan
sig = zeros(samplePoints,Ndetectors);
OCT_all = zeros(samplePoints/2+1,Ndetectors,Bscan_number);
for j_all = 1:Bscan_number
    for j = 1:Ndetectors
        disp(j)
        myname = ['sig1/',sig_type,'_',num2str(j),'_',num2str(j_all)];
        %% Load path lengths of detected photons DetS
        filename = sprintf('%s.bin',myname);
        if exist(filename,'file')~=0
            fid = fopen(filename, 'rb');
            Detsig = [fread(fid, 'double')];
            fclose(fid);
            sig_temp = Detsig(1:2:end)+1i*Detsig(2:2:end);
            sig(:,j) = sig_temp;
        end
    end
    sig = interp1(k_nolin,sig,k_lin);
    %% %%%%%%%
    %Adding noise
    
    noise1 = raylrnd(noise_amp/1.253*ones(size(sig))) .* exp(1i.*2.*pi.*rand(size(sig)));
    noise2 = raylrnd(noise_amp/1.253*ones(size(sig))) .* exp(1i.*2.*pi.*rand(size(sig)));
    ref_amp = 10*mean(mean(abs(sig),2));
    
    %Make the interference with the reference arm and calculate the intensity
    I = abs(sig + ref_amp+noise1).^2 - abs(sig - ref_amp+noise2).^2; % balance detection
    S_k =1;
    I = sqrt(S_k).*I;
    %% Processing the OCT signal
    k = (k_lin);
    %Apply low pass filter and hanning window
    ksampling = 2*pi/(k(1)-k(2));
    rawAline = I.*hanning(length(k));
    %Calculate Aline
    M10 = length(rawAline(:,1));
    OCT = abs(ifft(rawAline,M10));
    OCT = OCT(1:floor(length(OCT(:,1))/2)+1,:);
    OCT(2:end-1,:) = OCT(2:end-1,:);
    OCT_all(:,:,j_all) = OCT;
end
z = (0:M10/2)/M10.*ksampling;
x = linspace(-radius,radius,Ndetectors);
z = z/2;
%% OCTA
images = OCT_all;
X_size = samplePoints/2+1;
Y_size = Ndetectors;

oct_amplitude_temp = images(1:X_size,1:Y_size,:);
oct_amplitude = zeros(X_size+kernel_size-1,Y_size+kernel_size-1,Bscan_number);
for i=1:Bscan_number
    oct_amplitude(:,:,i) = padarray(oct_amplitude_temp(:,:,i),[(kernel_size-1)/2 (kernel_size-1)/2],'replicate','both');
end
oct_structure = mean(oct_amplitude,3);

bg_threshold = 0.0; %binary mask threshold (not used)
mask = ones(X_size+kernel_size-1,Y_size+kernel_size-1);
mask(oct_structure<bg_threshold)=0;

X_size1 = X_size+kernel_size-1;
Y_size1 = Y_size+kernel_size-1;
X_cm_map = X_size1 - kernel_size + 1;
Y_cm_map = Y_size1 - kernel_size + 1;
mask_used = mask(1:X_cm_map,1:Y_cm_map);

crop_intensity = oct_amplitude;
cm_map = zeros(X_cm_map,Y_cm_map,Bscan_number-1);

for k = 1:Bscan_number-1
    I_A = crop_intensity(:,:,k);
    I_B = crop_intensity(:,:,k+1);
    n = kernel_size.^2;
    for p = 1:X_cm_map
        for q = 1:Y_cm_map
            grid_A = I_A(p:p+kernel_size-1,q:q+kernel_size-1);
            grid_B = I_B(p:p+kernel_size-1,q:q+kernel_size-1);
            xy = grid_A.*grid_B;
            x2 = grid_A.^2;
            y2 = grid_B.^2;
            sum_xy = sum(xy,'all');
            sum_x = sum(grid_A,'all');
            sum_y = sum(grid_B,'all');
            r = (n.*sum_xy-sum_x.*sum_y)./(((n.*sum(x2,'all')-sum_x.^2).*(n.*sum(y2,'all')-sum_y.^2)).^0.5);
            cm_map(p,q,k) = r;
        end
    end
    
end
cor_map = cm_map;
cor_map_avg = mean(cor_map,3);
cor_map_avg = abs(cor_map_avg);
decor_map = 1-cor_map_avg;
decor_map_filtered = decor_map.*mask_used;

%% Display
mean1 = mean(oct_amplitude,3);
figure;
imagesc(x,z,10.*log10(mean1));
colormap gray
title('OCT B-scan (log scale)')
figure;
imagesc(x,z,decor_map_filtered);
colormap gray
title('OCTA')
%% Q-OCTA: print each vessel's velocity vs decorrelation
for id_vessel = 1:length(vessel_radius_all)
    d1 = decor_map_filtered;
    x1 = x*10;
    z1 = z*10/nm;
    x0 = vessel_center_xaxis_all(id_vessel)*10;% mm
    y0 = 0;
    z0 = vessel_center_zaxis_all(id_vessel)*10;% mm
    blood_radius = blood_radius_all(id_vessel)*10;% mm
    distances = zeros(length(z),length(x));
    for i= 1:length(z)
        for j=1:length(x)
            distances(i, j) = sqrt((z1(i) - z0)^2 + (x1(j) - x0)^2);
        end
    end
    z1 = z1(z1-z0<0);
    distances1 = distances(z1-z0<0,:);
    d1 = d1(z1-z0<0,:);
    RR = distances1(:)/blood_radius;
    d2 = d1(:);
    v1 = 2*u1*(1-RR.^2);
    v2 = v1(RR<=1);
    decor_map_filtered1 = d2(RR<=1);
    
    v21 = v2;
    decor_map_filtered2 = decor_map_filtered1;
    decor_map_filtered1(v21>velocity_threshold & decor_map_filtered2<decorrelation_threshold)=[];
    v2(v21>velocity_threshold & decor_map_filtered2<decorrelation_threshold)=[];
    figure
    data = [v2,decor_map_filtered1];
    
    x = data(:,1);
    y = data(:,2);
    
    edges = linspace(0, 1, numBins+1);
    binMeans = zeros(numBins, 1);
    binVars = zeros(numBins, 1);
    
    for i = 1:numBins
        binIndices = x > edges(i) & x <= edges(i+1);
        binData = y(binIndices);
        
        if ~isempty(binData)
            lowerBound = prctile(binData, lower_bound_num);
            upperBound = prctile(binData, upper_bound_num);
            
            filteredData = binData(binData >= lowerBound & binData <= upperBound);
            if ~isempty(filteredData)
                binMeans(i) = mean(filteredData);
                binVars(i) = std(filteredData);
            else
                binMeans(i) = NaN;
                binVars(i) = NaN;
            end
        else
            binMeans(i) = NaN;
            binVars(i) = NaN;
        end
    end
    
    binCenters = edges(1:end-1) + diff(edges)/2;
    binCenters = binCenters-binCenters(1);
    errorbar(binCenters, binMeans, binVars, '-o');
    xlabel('X Value');
    ylabel('Mean of Y');
    title('Binned Mean Y values with Variance Bars (25-75% Percentiles)');
    grid on;
end