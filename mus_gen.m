clear
close all
%% %%%%%%%%%%%%%%%%%%%%% Params settings DONE %%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dia = [4e-6,6e-6];% m, particle diameter for static medium and dynamic medium
pho = [0.00005,0.00025]; % number density for static medium and dynamic medium
samplePoints = 1024; % number of wavelength sampling
nm = 1.33; % backgroud refractive index
ns = 1.58; % particle refractive index
nang = 3601; % sampling angles number of SPF
lam_begin = 1250e-9; % m
lam_end = 1350e-9; % m
lam = linspace(lam_begin,lam_end,samplePoints); % wavelengths
%% %%%%%%%%%%%%%%%%%%%%% Params settings DONE %%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:2 % j=1: static; j=2: dynamic
    [~,ang,Miee,C,Q] = Mie(lam,dia(j),ns,nm,nang);
    Nb(j)=pho(j)*(1e18)*(4/3*pi*(dia(j)/2)^3); % concentration
    for i=1:length(lam)
        mus(i)=(C(i).sca)*pho(j)*1e16;
        r(i) = sqrt(Q(i).sca)*dia(j)/2*1e2;
    end
    save(['parameters/mus_',num2str(dia(j)*1e6),'_1300.mat'],'mus'); % cm, save mus
    save(['parameters/r_',num2str(dia(j)*1e6),'_1300.mat'],'r'); % cm, save r
end

function [S,ang,f,C1,Q1] = Mie(lam,dia,ns,nm,nang)
for i =1:length(lam)
    
    lambda = lam(i);    % vacuum wavelength
    conv = 1;           % convergence factor
    rad = dia/2.;           % sphere radius
    k = 2*pi/lambda*nm;    % wavenumber in medium n_m
    %% Calculate amplitude scattering matrix
    [S, C, ang,~] = calcmie(rad, ns, nm, lambda, nang, ...
        'ConvergenceFactor', conv);
    %% Calculate cross sections and efficiencies
    Q = getEfficiencies(C, rad(end), 3);
    C1(i) = C;
    Q1(i) = Q;
    %% test
    PP = (squeeze(abs(S(1,1,:).^2))+squeeze(abs(S(2,2,:).^2)))/pi/k/k/(Q.sca*rad(end)^2)/2;
    PP = PP/sum(PP);
    f(:,i) = PP;
end
end