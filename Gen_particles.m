clear
particle_dir = 'particles';
delete(fullfile(particle_dir, '*.bin')); % delete all old data in the new dictionary
static_particle_dir = [particle_dir,'/data_static_particles.bin'];
%% %%%%%%%%%%%%%%%%%%%%% Params settings SRART %%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('parameters/vessel_information.mat')
pho = [0.00005,0.00025]; %[Static, Dynamic],number density of particles, make sure it is the same as mus_gen.m
Nx = 400;% # of bins in each dimension of cube
Ny = 400;% # of bins in each dimension of cube
Nz = 200;% # of bins in each dimension of cube
binsize = 0.01; % mm,binsize, make sure it is the same as GenTissue.m
v=Nx*Ny*Nz;
u1 = 0.5; %mm/s, averaged velocity
D = 10; % mu_m^2/s, diffusion coefficient
delta_t = 0.01;% s, scanning rate
Bscan_number = 4; %All Bscan times, make sure it is the same as fullwaveMC_Bscan_OCTA.c
%% %%%%%%%%%%%%%%%%%%%%% Params settings DONE %%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Static particles: everywhere
X_range_all = Nx*binsize; % mm,X range
Y_range_all = Ny*binsize; % mm,Y range
Z_range_all = Nz*binsize; % mm,Z range
numPoints = int32(X_range_all*Y_range_all*Z_range_all*pho(1)*1e9);  % number of particles
points = zeros(numPoints, 4);
for i = 1:numPoints
    x = (2*rand-1) * X_range_all/2;
    y = (2*rand-1) * Y_range_all/2;
    z = rand * Z_range_all;
    ix = floor(Nx / 2 + x / binsize);
    iy = floor(Ny / 2 + y / binsize);
    iz = floor(z/binsize);
    i_all = iz * Ny * Nx + ix * Ny + iy;
    
    points(i, 1:3) = [x, y, z]/10;
    points(i, 4) = i_all;
end
points_sort = sortrows(points,4);
filename = static_particle_dir;
fileID = fopen(filename, 'wb');
fwrite(fileID, points_sort', 'double');
fclose(fileID);

%% Dynamic particles: only blood
for id_vessel = 1:length(vessel_radius_all)
    vessel_radius = vessel_radius_all(id_vessel)*10;% mm, make sure it is the same as GenTissue.m
    blood_radius = blood_radius_all(id_vessel)*10;% mm, blood_radius, make sure it is the same as GenTissue.m
    X_range = Nx*binsize; % mm,X range
    Y_range = Ny*binsize; % mm,Y range
    Z_range = vessel_radius*2; % mm,Z range
    x0 = vessel_center_xaxis_all(id_vessel)*10;% mm, make sure it is the same as GenTissue.m
    y0 = 0;
    z0 = vessel_center_zaxis_all(id_vessel)*10;% mm, make sure it is the same as GenTissue.m
    ymin = -Y_range/2;
    ymax = Y_range/2;
    xmin = x0-blood_radius;
    xmax = x0+blood_radius;
    zmin = z0-blood_radius;
    zmax = z0+blood_radius;
    %% Initial particles %%
    numPoints = int32(X_range*Y_range*Z_range*pho(2)*1e9);  % number of particles
    points = zeros(numPoints, 4);
    for i = 1:numPoints
        x = (2*rand-1) * X_range/2;
        y = (2*rand-1) * Y_range/2;
        z = rand * Z_range+z0-vessel_radius;
        ix = floor(Nx / 2 + x / binsize);
        iy = floor(Ny / 2 + y / binsize);
        iz = floor(z/binsize);
        i_all = iz * Ny * Nx + ix * Ny + iy;
        
        points(i, 1:3) = [x, y, z]/10;
        points(i, 4) = i_all;
    end
    points_sort = sortrows(points,4);
    
    filename =[particle_dir,'/data_dynamic_particles_init_vessel',num2str(id_vessel),'.bin'];
    fileID = fopen(filename, 'wb');
    fwrite(fileID, points_sort', 'double');
    fclose(fileID);
    
    
    %% Dynamic flow %%
    filename = [particle_dir,'/data_dynamic_particles_init_vessel',num2str(id_vessel),'.bin'];
    fileID = fopen(filename);
    A = fread(fileID,'double');
    fclose(fileID);
    A1 = reshape(A,4,[])';
    A2 = A1;
    A1(:,1:3) = A1(:,1:3)*10;
    R = sqrt((A1(:,1)-x0).^2 + (A1(:,3)-z0).^2)/blood_radius;
    RRR = R(R<=1);
    filename = [particle_dir,'/data_dynamic_particles_init_vessel',num2str(id_vessel),'.bin'];
    fileID = fopen(filename, 'wb');
    fwrite(fileID, A2(R<=1,:)', 'double');
    fclose(fileID);
    for Bscan_times = 1:Bscan_number-1
        if Bscan_times==1
            filename = [particle_dir,'/data_dynamic_particles_init_vessel',num2str(id_vessel),'.bin'];
        else
            filename = [particle_dir,'/data_dynamic_particles_vessel',num2str(id_vessel),'_',num2str(Bscan_times),'.bin'];
        end
        fileID = fopen(filename);
        A = fread(fileID,'double');
        fclose(fileID);
        A1 = reshape(A,4,[])';
        A1(:,1:3) = A1(:,1:3)*10;
        A2 = A1;
        R = sqrt((A1(:,1)-x0).^2 + (A1(:,3)-z0).^2)/blood_radius;
        u = 2*u1*(1-R.^2);
        rnd_br = sqrt(2*D*delta_t*1)*1e-3*randn(size(A1,1),3);
        A2(R<=1,2) = A2(R<=1,2)+u(R<=1)*delta_t*1;
        A2(R<=1,1:3) = A2(R<=1,1:3)+rnd_br(R<=1,:);
        points = zeros(length(A2), 4);
        for i = 1:length(A2)
            x = A2(i,1);
            y = A2(i,2);
            z = A2(i,3);
            if x>xmax
                x=xmax;
            end
            if y>ymax
                y=ymax;
            end
            if z>zmax
                z=zmax;
            end
            if x<xmin
                x=xmin;
            end
            if y<ymin
                y=ymin;
            end
            if z<zmin
                z=zmin;
            end
            ix = floor(Nx / 2 + x / binsize);
            iy = floor(Ny / 2 + y / binsize);
            iz = floor(z/binsize);
            i_all = iz * Ny * Nx + ix * Ny + iy;
            points(i, 1:3) = [x, y, z]/10;
            points(i, 4) = i_all;
        end
        R1 = sqrt((points(:,1)*10-x0).^2 + (points(:,3)*10-z0).^2)/blood_radius;
        points2 = points(R1<=1,:);
        points_sort = sortrows(points2,4);
        
        filename = [particle_dir,'/data_dynamic_particles_vessel',num2str(id_vessel),'_',num2str(Bscan_times+1),'.bin'];
        fileID = fopen(filename, 'wb');
        fwrite(fileID, points_sort', 'double');
        fclose(fileID);
    end
end

for i = 1:Bscan_number
    data = [];
    if i == 1
        filename_towrite = [particle_dir,'/data_dynamic_particles_init.bin'];
        
    else
        filename_towrite = [particle_dir,'/data_dynamic_particles_',num2str(i),'.bin'];
    end
    for j = 1:length(vessel_radius_all)
        if i == 1
            filename_toread = [particle_dir,'/data_dynamic_particles_init_vessel',num2str(j),'.bin'];
        else
            filename_toread = [particle_dir,'/data_dynamic_particles_vessel',num2str(j),'_',num2str(i),'.bin'];
        end
        fileID = fopen(filename_toread);
        data_temp = fread(fileID,'double');
        
         data_temp = reshape(data_temp,4,[])';
         data = [data;data_temp];
        fclose(fileID);
    end
    points_sort = sortrows(data,4);
    fileID = fopen(filename_towrite, 'ab');
    fwrite(fileID, points_sort', 'double');
    fclose(fileID);
end
