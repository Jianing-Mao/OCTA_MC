clear
close all
dia = 4e-6; % particle diameter for static medium
% dia = 6e-6; % particle diameter for dynamic medium
samplePoints = 1024; % number of wavelength sampling
nm=1.33; % backgroud refractive index
ns = 1.58; % particle refractive index
nang = 3601; % sampling angles number of SPF 
lam_begin = 1250e-9;
lam_end = 1350e-9;
lam = linspace(lam_begin,lam_end,samplePoints); % wavelengths
pho = 0.00005; % number density for static medium
% pho = 0.00025; % number density for dynamic medium
[~,ang,Miee,C,Q] = Mie(lam,dia,ns,nm,nang);
k0 = 2*pi/(1260/1e9)*nm;
Nb=pho*(1e18)*(4/3*pi*(dia/2)^3); % concentration
for i=1:length(lam)
us(i)=(C(i).sca)*pho*1e16;
r(i) = sqrt(Q(i).sca)*dia/2*1e2;
end
eval(['us_', num2str(dia*1e6), '=', 'us', ';']);
save(['us_',num2str(dia*1e6),'_1300.mat'],['us_', num2str(dia*1e6)]); %save
eval(['r_', num2str(dia*1e6), '=', 'r', ';']);
save(['r_',num2str(dia*1e6),'_1300.mat'],['r_', num2str(dia*1e6)]); %save

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