%% Create the files for the simulation
clear
format compact
clc
home
%% load the desired mus and TT parameters of different particle diameters
load us_4_1300.mat % scattering coefficient of 4 mu_m diameter particles
load us_6_1300.mat % scattering coefficient of 6 mu_m diameter particles
load r_4_1300.mat % efficient particle radius of 6 mu_m diameter particles
load r_6_1300.mat % efficient particle radius of 6 mu_m diameter particles
para4 = load('parameters_4.txt'); % TT SPF parameters of 4 mu_m diameter particles
para6 = load('parameters_6.txt'); % TT SPF parameters of 6 mu_m diameter particles

%% set paras
%%% USER CHOICES %%%%%%%%
SAVEON      = 1;        % 1 = save myname_T.bin, myname_H.mci
samplePoints = 1024; %number of wavelength sampling
myname      = 'infi';% name for files: myname_T.bin, myname_H.mci
Nphotons    = 15000;      	% number of photons
Nx          = 400;    	% # of bins in each dimension of cube
Ny          = 400;    	% # of bins in each dimension of cube
Nz          = 200;    	% # of bins in each dimension of cube
binsize     = 0.001;     	% size of each bin, eg. [cm]

n = 1.33; % refractive index

% Set detection parameter
radius      = 0.2;     % Half width of the BScan
Ndetectors  = 512;      % Number of Aline per BScan, 1 for Aline simulation
beamw = 0.2;
flens = 3.6;
% Cos of the accepted angle
cos_accept  = cos(atan(beamw/2/flens)); % Cos of the accepted angle
z_focus = 0.08;% Focal depth

% Bias scattering parameter
p  = 0.5;
dx = binsize;
dy = binsize;
dz = binsize;
%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify Monte Carlo parameters
x  = ([0:Nx]'-Nx/2)*dx;
y  = ([0:Ny]'-Ny/2)*dy;
z  = [0:Nz]'*dz;
zmin = min(z);
zmax = max(z);
xmin = min(x);
xmax = max(x);
%%%%%%%%%%%%
%% Create Sample
%%% number of mediums
Nt = 9;
%%% set your structure and parameters for each wavelength
lam_begin = 1250e-7;
lam_end = 1350e-7;
lam = linspace(lam_begin,lam_end,samplePoints); % wavelengths
start_depth = 0.01;% cm
vessel_center_depth = 0.08;% cm
vessel_radius = 0.07;% cm
blood_radius = 0.03;% cm
Z_range = vessel_radius*2;% cm
Y_range = Ny*binsize;
X_range = Nx*binsize;

%Please use medium number=5 as the dynamic blood;medium number=4 as the vessel wall;
for jj=1:samplePoints
    myname = ['infi',num2str(jj)]
    waist  = (lam(jj)/n)/pi/(beamw/2/flens);    % Width of the beam at the imaging lens
    zr = pi*waist^2/(lam(jj)/n);
    det_radius = waist*sqrt(1+(z_focus/zr)^2);
    for i=1:Nt
        if i == 1
            nrv(i)   = n;
        elseif i == 4
            nrv(i) = n; % refraction index (set uniform for every medium)
            muav(i) = 0; % absorption coefficient
            musv(i)= us_4(jj); % scattering coefficient
            gv(i) = 0.9; % anisotropy
            gf(i) = para4(jj,1); % forward anisotropy
            gb(i) = para4(jj,2); % backward anisotropy
            alf(i) = para4(jj,3); % forward enhancing
            alb(i) = para4(jj,4); % backward enhancing
            C(i) = para4(jj,5); % balance factor
            r(i) = r_4(jj); %efficient r
        elseif i==5
            nrv(i) = n;
            muav(i) = 0;
            musv(i)= us_6(jj);
            gv(i) = 0.9;
            gf(i) = para6(jj,1);
            gb(i) = para6(jj,2);
            alf(i) = para6(jj,3);
            alb(i) = para6(jj,4);
            C(i) = para6(jj,5);
            r(i) = r_6(jj); %efficient r
        else
            nrv(i) = 1.4+rand;
            muav(i) = 0;
            musv(i)= 0;
            gv(i) = 0;
            gf(i) = 0;
            gb(i) = 0;
            alf(i) = 0;
            alb(i) = 0;
            C(i) = 0;
            r(i) = 0;
        end
    end
    T = double(zeros(Ny,Nx,Nz));
    
    T = T + 1;      %
    
    zsurf = 0.0000;  % position of surface
    %%% set the structure
    for iz=1:Nz % for every depth z(iz)
           if iz<=round(start_depth/dz)
            T(:,:,iz) =1;
           end
        
    xc      = 0;            % [cm], center of blood vessel
    zc      = vessel_center_depth;     	% [cm], center of blood vessel
    wallradius  = vessel_radius;      	% blood vessel radius [cm]
    for ix=1:Nx
            xd = x(ix) - xc;	% vessel, x distance from vessel center
            zd = z(iz) - zc;   	% vessel, z distance from vessel center                
            rr  = sqrt(xd^2 + zd^2);	% r from vessel center
            if (rr<=wallradius)     	% if r is within vessel
                T(:,ix,iz) =4; % vessel wall
            end

    end %ix    
    xc      = 0;            % [cm], center of blood vessel
    zc      = vessel_center_depth;     	% [cm], center of blood vessel
    vesselradius  = blood_radius;      	% blood vessel radius [cm]
    for ix=1:Nx
            xd = x(ix) - xc;	% vessel, x distance from vessel center
            zd = z(iz) - zc;   	% vessel, z distance from vessel center                
            rr  = sqrt(xd^2 + zd^2);	% r from vessel center
            if (rr<=vesselradius)     	% if r is within vessel
                T(:,ix,iz) =5; % blood
            end

    end %ix
    end % iz
    
    %%
    if SAVEON
        tic
        % convert T to linear array of integer values, v(i)i = 0;
        if jj==1
        v = uint8(reshape(T,Ny*Nx*Nz,1));
        end
        %% WRITE FILES
        % Write myname_H.mci file
        %   which contains the Monte Carlo simulation parameters
        %   and specifies the tissue optical properties for each tissue type.
        commandwindow
        disp(sprintf('--------create %s --------',myname))
        filename = sprintf('%s_H.mci',myname);
        fid = fopen(filename,'w');
        % run parameters
        fprintf(fid,'%0.2f\n',Nphotons);
        fprintf(fid,'%0.4f\n',p);
        fprintf(fid,'%0.4f\n',Ndetectors);
        fprintf(fid,'%0.6f\n',det_radius);
        fprintf(fid,'%0.6f\n',cos_accept);
        fprintf(fid,'%d\n'   ,Nx);
        fprintf(fid,'%d\n'   ,Ny);
        fprintf(fid,'%d\n'   ,Nz);
        fprintf(fid,'%0.4f\n',dx);
        fprintf(fid,'%0.4f\n',dy);
        fprintf(fid,'%0.4f\n',dz);
        % launch parameters
        fprintf(fid,'%0.4f\n',radius);
        fprintf(fid,'%0.4f\n',zsurf);
        % tissue optical properties
        fprintf(fid,'%d\n',Nt);
        for i=1:Nt
            fprintf(fid,'%0.6f\n',muav(i));
            fprintf(fid,'%0.6f\n',musv(i));
            fprintf(fid,'%0.6f\n',gv(i));
            fprintf(fid,'%0.6f\n',nrv(i));
            fprintf(fid,'%0.6f\n',gf(i));
            fprintf(fid,'%0.6f\n',gb(i));
            fprintf(fid,'%0.6f\n',alf(i));
            fprintf(fid,'%0.6f\n',alb(i));
            fprintf(fid,'%0.6f\n',C(i));
            fprintf(fid,'%0.6f\n',r(i));
        end
        fprintf(fid,'%0.4f\n',z_focus);
        fprintf(fid,'%0.4f\n',waist);
        fprintf(fid,'%0.4f\n',zr);
        fclose(fid);
        
        %% write myname_T.bin file
        if jj==1
            filename = sprintf('%s_T.bin',myname);
            disp(['create ' filename])
            fid = fopen(filename,'wb');
            fwrite(fid,v,'uint8');
            fclose(fid);
        end
        toc
    end % SAVEON
end
%% Look at structure of Tzx at iy=Ny/2
Txzy = shiftdim(T,1);   % Tyxz --> Txzy
Tzx  = Txzy(:,:,Ny/2)'; % Tzx

%%
figure('color', 'white'); clf
sz = 12;  fz = 10;
imagesc(x,z,Tzx,[1 Nt])
hold on
set(gca,'fontsize',sz)
xlabel('x [cm]')
ylabel('z [cm]')
title('\rm Tissue')
colorbar
cmap = makecmap(Nt);
colormap(cmap)
set(colorbar,'fontsize',1)
% label colorbar
zdiff = zmax-zmin;
%%%
text(xmax,zmin - zdiff*0.06, 'Tissue types','fontsize',fz)
axis equal image
axis([xmin xmax zmin zmax])

savefig(strcat(myname, '.fig'));
disp('done')


