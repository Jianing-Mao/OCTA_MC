%% params
clear
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%% Params settings START %%%%%%%%%%%%%%%%%%%%%%%%%%%
dia = [4e-6,6e-6];% m, particle diameter for static medium and dynamic medium
eps = 1e-9;
samplePoints = 1024; % number of wavelength sampling
nm = 1.33; % backgroud refractive index
ns = 1.58; % particle refractive index
nang = 3601; % sampling angles number of SPF 
lam_begin = 1250e-9;
lam_end = 1350e-9;
lam = linspace(lam_begin,lam_end,samplePoints); % wavelengths
%% %%%%%%%%%%%%%%%%%%%%% Params settings DONE %%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:2 % j=1: static; j=2: dynamic
[~,ang,Miee,C] = Mie(lam,dia(j),ns,nm,nang);
ang_pi = ang/180*pi;
x_final = zeros(samplePoints,5); % parameters of TT
atten = [linspace(1,1,(nang-1)/2),linspace(1,5,(nang-1)/2+1)]; % scaling
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective','MaxIterations',100000);
%% Fit
for i=1:samplePoints
    disp(i)
    ydata = log10(Miee(:,i)').*atten;
    xdata = ang_pi;
    fun = @(x,xdata)log10(TT(xdata,x(1),x(2),x(3),x(4),x(5))).*atten;
    x0 = [0.5,0.5,1,1,0.99]; %[gf,gb,alpha_f,alpha_b,C]
    lb = [0+eps,0+eps,-0.5+eps,-0.5+eps,0];
    ub = [1,1,10,10,1];
    x = lsqcurvefit(fun,x0,xdata,ydata,lb,ub,options);
    x_final(i,:) = x;
end
dlmwrite(['parameters/parameters_',num2str(dia(j)*1e6),'.txt'],x_final,'precision',20)
end

%% function
function [S,ang,f,C1] = Mie(lam,dia,ns,nm,nang)
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
    %% test
    PP = (squeeze(abs(S(1,1,:).^2))+squeeze(abs(S(2,2,:).^2)))/pi/k/k/(Q.sca*rad(end)^2)/2;
    f(:,i) = PP/max(PP);
end
end

function f = TT(ang,gf,gb,alpha_f,alpha_b,C)
Kf = 1 / pi * alpha_f * gf * (1 - gf^2)^(2 * alpha_f) / ((1 + gf)^(2 * alpha_f) - (1 - gf)^(2 * alpha_f));
Kb = 1 / pi * alpha_b * gb * (1 - gb^2)^(2 * alpha_b) / ((1 + gb)^(2 * alpha_b) - (1 - gb)^(2 * alpha_b));
ang_re = flip(ang);
ff1 = (1 + (gf * gf) - 2 * gf * cos(ang));
ff2 = ff1.^( -(alpha_f + 1));
ff = Kf * ff2;
fb1 = (1 + gb*gb - 2 * gb * cos(ang_re));
fb2 = fb1.^(-(alpha_b + 1));
fb = Kb * fb2;
f = C * ff + (1 - C) * fb;
f = f/max(f);
end