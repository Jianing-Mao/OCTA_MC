clear
%% Parameters %%
Nx = 400;% # of bins in each dimension of cube
Ny = 400;% # of bins in each dimension of cube
Nz = 200;% # of bins in each dimension of cube
binsize = 0.01; % mm,binsize, make sure it is the same as GenTissue.m
vessel_radius = 0.7;% mm, make sure it is the same as GenTissue.m
blood_radius = 0.3;% mm, blood_radius, make sure it is the same as GenTissue.m
X_range = Nx*binsize; % mm,X range, make sure it is the same as GenTissue.m
Y_range = Ny*binsize; % mm,Y range, make sure it is the same as GenTissue.m
Z_range = vessel_radius*2; % mm,Z range, make sure it is the same as GenTissue.m
start_depth = 0.1;% mm,start_depth, make sure it is the same as GenTissue.m
v=Nx*Ny*Nz;
u1 = 0.5; %mm/s, averaged velocity
D = 10; % mu_m^2/s, diffusion coefficient
delta_t = 0.01;% s, scanning rate
x0 = 0;
y0 = 0;
z0 = 0.8;% mm, vessel_center_depth, make sure it is the same as GenTissue.m
ymin = -Y_range/2;
ymax = Y_range/2;
xmin = x0-blood_radius;
xmax = x0+blood_radius;
zmin = z0-blood_radius;
zmax = z0+blood_radius;
Bscan_number = 4; %All Bscan times, make sure it is the same as fullwaveMC_Bscan_OCTA.c
pho = [0.00005,0.00025]; %[Static, Dynamic],number density of particles, make sure it is the same as mus_gen.m

%% Initial particles %%
for type = 1:2
    numPoints = X_range*Y_range*Z_range*pho(type)*1e9;  % number of particles
    points = zeros(numPoints, 4);  
    for i = 1:numPoints
     
        x = (2*rand-1) * X_range/2;
        y = (2*rand-1) * Y_range/2;
        z = rand * Z_range+start_depth;
        ix = floor(Nx / 2 + x / binsize);
        iy = floor(Ny / 2 + y / binsize);
        iz = floor(z/binsize);
        i_all = iz * Ny * Nx + ix * Ny + iy;
   
        points(i, 1:3) = [x, y, z]/10;
        points(i, 4) = i_all;
    end
    points_sort = sortrows(points,4);
    if type==1
        filename = 'data_static_particles.bin';
        fileID = fopen(filename, 'wb');
        fwrite(fileID, points_sort', 'double');
        fclose(fileID);
    else
        filename = 'data_dynamic_particles_init.bin';
        fileID = fopen(filename, 'wb');
        fwrite(fileID, points_sort', 'double');
        fclose(fileID);
    end
end

%% Dynamic flow %%
filename = 'data_dynamic_particles_init.bin';
fileID = fopen(filename);
A = fread(fileID,'double');
A1 = reshape(A,4,[])';
A2 = A1;
A1(:,1:3) = A1(:,1:3)*10;
R = sqrt((A1(:,1)-x0).^2 + (A1(:,3)-z0).^2)/blood_radius;
RRR = R(R<=1);
filename = 'data_dynamic_particles_init.bin';
fileID = fopen(filename, 'wb');
fwrite(fileID, A2(R<=1,:)', 'double');
fclose(fileID);
for Bscan_times = 1:Bscan_number-1
    disp(Bscan_times)
    if Bscan_times==1
        filename = 'data_dynamic_particles_init.bin';
    else
        filename = ['data_dynamic_particles_',num2str(Bscan_times),'.bin'];
    end
fileID = fopen(filename);
A = fread(fileID,'double');
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

filename = ['data_dynamic_particles_',num2str(Bscan_times+1),'.bin'];
fileID = fopen(filename, 'wb');
fwrite(fileID, points_sort', 'double');
fclose(fileID);
end

