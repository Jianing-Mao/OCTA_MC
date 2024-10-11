#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <complex.h>
#include <omp.h>

#define Ntiss 19   /* Number of tissue types. */
#define STRLEN 32  /* String length. */
#define ls 1.0E-10 /* Moving photon a little bit off the voxel face */
#define PI 3.1415926
#define LIGHTSPEED 2.997925E10 /* in vacuo speed of light [cm/s] */
#define ALIVE 1                /* if photon not yet terminated */
#define DEAD 0                 /* if photon is to be terminated */
#define THRESHOLD 0.01         /* used in roulette */
#define CHANCE 0.1             /* used in roulette */
#define Boolean char
#define SQR(x) (x * x)
#define SIGN(x) ((x) >= 0 ? 1 : -1)
#define RandomNum (double)RandomGen(1, 0, NULL) /* Calls for a random number. */
#define COS90D 1.0E-6                           /* If cos(theta) <= COS90D, theta >= PI/2 - 1e-6 rad. */
#define ONE_MINUS_COSZERO 1.0E-12
#define Sep_photons 1000
#define ORDER 4
#define nang 36001
#define nang_is 18001

#define samplePoints 1024
#define lambda_start 1250e-7
#define lambda_end 1350e-7

typedef struct {
    double x, y, z;
} Point3D;

typedef struct {
    Point3D start;
    Point3D end;
} Line3D;
/* DECLARE FUNCTIONS */
double RandomGen(char Type, long Seed, long *Status);
/* Random number generator */
Boolean SameVoxel(double x1, double y1, double z1, double x2, double y2, double z2,
                  double dx, double dy, double dz);
/* Asks,"In the same voxel?" */
double max2(double a, double b);
double min2(double a, double b);
double min3(double a, double b, double c);
double FindVoxelFace2(int *fdir, double x1, double y1, double z1, int *det_num, int Pick_det,
                      double detx, double det_radius, double det_z, double cos_accept,
                      int Ndetectors, double dx, double dy, double dz, double ux, double uy, double uz);
double TT_final1(double gf, double gb, double alpha_f, double alpha_b, double C, double ang);
double TT_is_final1(double gf, double gb, double alpha_f, double alpha_b, double C, double ang);
void fit_polynomial(double *x, double *y, int num_points, double *coeffs);
double polyval(double *coeffs, int order, double x);
void linspace(double start, double stop, int num, double *arr);
void linear_interpolate(double *x, double *y, int n, double *x_interp, double *y_interp, int m);
void write_complex_array_to_binfile(char *filename, double complex *data, int size);
double distanceToLineSegment(Point3D point, Point3D start, Point3D end);
void reshapeArray(const double* input, int rows, int cols, double*** output);
void findIndices(const long* array, int size, long target, int* lowerIndex, int* upperIndex);
double RFresnel(double n1,double n2,double ca1,double *ca2_Ptr);
void CrossUpOrNot(double ni, double nt);
void CrossDnOrNot(double ni, double nt);
void CrossOutOrNot(double ni, double nt);
void CrossInOrNot(double ni, double nt);
void CrossLeftOrNot(double ni, double nt);
void CrossRightOrNot(double ni, double nt);
void generateRandomGaussian(double *x, double *y,double r1, double r2);
void generateGaussianDirection(double x, double y, double z, double *ux, double *uy, double *uz, double zr, double zfocus, double detx);
void calculateReflection(double* reflectedX, double* reflectedY, double* reflectedZ, 
                         double incidentX, double incidentY, double incidentZ, 
                         double pointAX, double pointAY, double pointAZ, 
                         double pointBX, double pointBY, double pointBZ);
double cos_is_final1(double gf, double gb, double alpha_f, double alpha_b, double C);
double cos_final1(double gf, double gb, double alpha_f, double alpha_b, double C);


/* Propagation parameters */
double x, y, z,x_last, y_last, z_last;       /* photon position */
double ux, uy, uz;    /* photon trajectory as cosines */
double uxx, uyy, uzz; /* temporary values used during SPIN */
double s;             /* step sizes. s = -log(RND)/mus [cm] */
double sleft;         /* dimensionless */
double costheta;      /* cos(theta) */
double sintheta;      /* sin(theta) */
double cospsi;        /* cos(psi) */
double sinpsi;        /* sin(psi) */
double psi;           /* azimuthal angle */
long i_photon;        /* current photon */
double W;             /* photon weight */
double absorb;        /* weighted deposited in a step due to absorption */
short photon_status;  /* flag = ALIVE=1 or DEAD=0 */
Boolean sv;           /* Are they in the same voxel? */
int fdir;             /* 1:x,2:y,3:z */

/* other variables */
double mua;      /* absorption coefficient [cm^-1] */
double mus;      /* scattering coefficient [cm^-1] */
double g;        /* anisotropy [-] */
double nr;       /* refractive index [-] RMT */
double Nphotons; /* number of photons in simulation */
double effr;

double gf;      /* forward anisotropy [-] */
double gb;      /* backward anisotropy [-] */
double alpha_f; /* forward enhance [-] */
double alpha_b; /* backward enhance [-] */
double C;       /* balance [-] */

/* launch parameters */
float ux0, uy0, uz0;
float radius;
float zfocus,zr,waist,zfocus_used;

/* dummy variables */
double rnd;                 /* assigned random value 0-1 */
double rr, phi;             /* dummy values */
long i, j, NN, Nyx,i_particle_center,i_particle;         /* dummy indices */
double tempx, tempy, tempz; /* temporary variables, used during photon step. */
int ix, iy, iz;             /* Added. Used to track photons */
double temp;                /* dummy variable */
int ix_particle, iy_particle, iz_particle;

int surfflag; /* surface flag: 0 = photon inside tissue, 1 = escaped outside tissue */

/* mcxyz bin variables */
float dx, dy, dz;   /* bin size [cm] */
int Nx, Ny, Nz, Nt; /* # of bins */
float xs, ys, zs;   /* launch position */
float zsurf, Rd;

/* time */
float time_min; // Requested time duration of computation.
time_t now;
double start_time, finish_time, temp_time, start_time_all, end_time_all; /* for clock() */

/* tissue parameters */
char tissuename[50][32];
float muav[Ntiss]; // muav[0:Ntiss-1], absorption coefficient of ith tissue type
float musv[Ntiss]; // scattering coeff.
float nrv[Ntiss];  // refractive index
float gv[Ntiss];   // anisotropy of scattering
float gf_tt[Ntiss];
float gb_tt[Ntiss];
float alf_tt[Ntiss];
float alb_tt[Ntiss];
float C_tt[Ntiss];
float TT[Ntiss][nang] = {{0}};
float TT_is[Ntiss][nang_is] = {{0}};
float max_TT_v[Ntiss] = {0};
float max_TT_is_v[Ntiss] = {0};
float eff_r[Ntiss];

int det_num;             // photon not detected yet/
int first_bias_done;  // photon not biased back - scattered yet
int cont_exist;       // no split generated yet // check if a continuing photon packet exists
int haven_scatter;
double L_current;        // photon 's initial likelihood
double s_total,s_total_sc;          // photon 's initial path length
double z_max;            // photon 's initial depth reached
int Ndetectors;          // Number of source/detector pairs to represent conducting A-scans
int Pick_det;            // index of randomly picked source/detector pair
double detx, dety, detz; // position of detector
double det_radius;       // radius of detector
double cos_accept;       // acceptance angle
double a_coef;           // Parameter of the bias function
double costheta_S;
double costheta_B;
double sintheta_B;
double vx, vy, vz;
double upx, upy, upz;
double L_cont;
long i_cont;
double W_cont;
double tempslen_cont;
double x_cont, y_cont, z_cont, x_old_cont, y_old_cont, z_old_cont;
double ux_cont, uy_cont, uz_cont;
double s_total_cont, num_s_cont;
double z_max_cont;
double p; // parameter of chance of biased forward-scattering
double det_z;
double f_TT, f_B;
long c_photon,c_photon_ss,c_photon_ms; // count collected photons
int *Detsnum = NULL;
float *DetW = NULL, *DetL = NULL, *DetS = NULL,*DetW_ss = NULL, *DetL_ss = NULL, *DetS_ss = NULL,*DetW_ms = NULL, *DetL_ms = NULL, *DetS_ms = NULL;
// float *DetW = NULL, *DetL = NULL, *DetS = NULL;
double temp11, temp22, L_temp;
int k = 0;
int nnnn = 0;
int type, split_num, split_num_cont,type_old,type_old_cont;
char *v = NULL;
double *particle_temp = NULL;
int num_s;

int ii = 0;
int back_pick = 0;
int find_sample = 0;
double max_TT = 0, ran1 = 0, ran2 = 0, Sample = 0, max_TT_is = 0;
int iindex = 0;
int scatter_find = 0;
double scatter_dis = 0;
double MM[samplePoints][samplePoints] = {{0}};
double fit[samplePoints][samplePoints] = {{0}};
double MM_ss[samplePoints][samplePoints] = {{0}};
double MM_ms[samplePoints][samplePoints] = {{0}};
double fit_BDR[samplePoints][samplePoints] = {{0}};


int temp_lin, temp_lin2, temp_lin3;
double lambda[samplePoints], k_nolin[samplePoints], k_lin[samplePoints];
double temp_x[samplePoints];
double coeffs[ORDER + 1];
double fit_temp[samplePoints], sum_MM, sum_fit, fit_temp_MM[samplePoints];
double fit_trans[samplePoints][samplePoints] = {{0}};

int c_photon_mask[samplePoints] = {0};
int c_photon_mask_ss[samplePoints] = {0};
int c_photon_mask_ms[samplePoints] = {0};
int particle_finde_index=-1; int particle_finde_index_center = -1;



int main(int argc, const char *argv[])
{
    //    printf("argc = %d\n",argc);
    if (argc == 0)
    {
        printf("which will load the files name_H.mci and name_T.bin\n");
        printf("and run the Monte Carlo program.\n");
        return 0;
    }
    
    /* Input/Output */
    char myname[STRLEN];
    // char Miename[STRLEN];
    char Miename[STRLEN],Miename_ss[STRLEN],Miename_ms[STRLEN];
    char Miename2[STRLEN];
    // Holds the user's choice of myname, used in input and output files.
    char filename[STRLEN]; // temporary filename for writing output.
    FILE *fid = NULL;      // file ID pointer
    int B_scan_num = atoi(argv[2]);


////////////////////////////////////////////
/*Load particles in the static medium*/
    char particle[STRLEN];
    strcpy(particle, "particles/data_static_particles.bin");
    fid = fopen(particle,"rb");
    fseek(fid, 0, SEEK_END);
    long fileSize = ftell(fid);
    fseek(fid, 0, SEEK_SET);
    long numDoubles = fileSize / sizeof(double);
    particle_temp = (double *)malloc(fileSize);
    fread(particle_temp, sizeof(double), numDoubles , fid);
    double** particles;
    int particle_num = numDoubles/4;
    reshapeArray(particle_temp, particle_num, 4, &particles);
    fclose(fid);
    free(particle_temp);
    printf("%f\n",particles[0][0]);
    printf("%f\n",particles[0][1]);
    printf("%f\n",particles[0][2]);
    printf("%f\n",particles[0][3]);
    long* particle_index;
    particle_index = (long*)malloc(particle_num * sizeof(long));
    for(int temp_i = 0;temp_i<particle_num;temp_i++){
        particle_index[temp_i] = (long)(particles[temp_i][3]);
    }
      

////////////////////////////////////////////
for(int iter_whole=1;iter_whole<=B_scan_num;iter_whole++){
    double complex sig[samplePoints] = {0 + 0 * I}; // Signals of all-photon settings
    double complex sig_ss[samplePoints] = {0 + 0 * I}; // Signals of single scattered-photon settings
    double complex sig_ms[samplePoints] = {0 + 0 * I}; // Signals of multiple scattered-photon settings
    char num_iterall[20];
    sprintf(num_iterall, "%d", iter_whole);

/*Load particles in the dynamic medium*/
    char particle_center[STRLEN];
    if(iter_whole==1)
    strcpy(particle_center, "particles/data_dynamic_particles_init.bin");
    else{
        strcpy(particle_center, "particles/data_dynamic_particles_");
        strcat(particle_center, num_iterall);
        strcat(particle_center, ".bin");
    }
    fid = fopen(particle_center,"rb");
    fseek(fid, 0, SEEK_END);
    long fileSize = ftell(fid);
    fseek(fid, 0, SEEK_SET);
    long numDoubles = fileSize / sizeof(double);
    particle_temp = (double *)malloc(fileSize);
    fread(particle_temp, sizeof(double), numDoubles , fid);
    double** particles_center;
    int particle_center_num = numDoubles/4;
    reshapeArray(particle_temp, particle_center_num, 4, &particles_center);
    fclose(fid);
    free(particle_temp);
    printf("%f\n",particles_center[0][0]);
    printf("%f\n",particles_center[0][1]);
    printf("%f\n",particles_center[0][2]);
    printf("%f\n",particles_center[0][3]);
    long* particle_center_index;
    particle_center_index = (long*)malloc(particle_center_num * sizeof(long));
    for(int temp_i = 0;temp_i<particle_center_num;temp_i++){
        particle_center_index[temp_i] = (long)(particles_center[temp_i][3]);
    }

    
    char buf[32];          // buffer for reading header.dat
    int p_num = atoi(argv[1]);
    int rndn = atoi(argv[1]);

    strcpy(Miename, "sig1/sig_"); // acquire name from argument of function call by user
     strcat(Miename, argv[1]);
    strcat(Miename, "_");
     strcat(Miename, num_iterall);
    strcat(Miename, ".bin");

    strcpy(Miename_ss, "sig1/sig_ss_"); // acquire name from argument of function call by user
     strcat(Miename_ss, argv[1]);
     strcat(Miename_ss, "_");
     strcat(Miename_ss, num_iterall);
    strcat(Miename_ss, ".bin");

    strcpy(Miename_ms, "sig1/sig_ms_"); // acquire name from argument of function call by user
     strcat(Miename_ms, argv[1]);
     strcat(Miename_ms, "_");
     strcat(Miename_ms, num_iterall);
    strcat(Miename_ms, ".bin");

    double z_min;
    DetS = malloc(sizeof(float)); // photon path length
    DetL = malloc(sizeof(float)); // likelihood ratio
    DetW = malloc(sizeof(float)); // photon weight

    DetS_ss = malloc(sizeof(float)); // photon path length
    DetL_ss = malloc(sizeof(float)); // likelihood ratio
    DetW_ss = malloc(sizeof(float)); // photon weight

    DetS_ms = malloc(sizeof(float)); // photon path length
    DetL_ms = malloc(sizeof(float)); // likelihood ratio
    DetW_ms = malloc(sizeof(float)); // photon weight

    start_time_all = clock();
    for (int aa = 1; aa <= samplePoints; aa++)
    {
        strcpy(myname, "settings/infi"); // acquire name from argument of function call by user
        char num[20];
        sprintf(num, "%d", aa);
        strcat(myname, num);
        printf("name = %s\n", myname);

        /**** INPUT FILES *****/
        /* IMPORT myname_H.mci */
        strcpy(filename, myname);
        strcat(filename, "_H.mci");
        fid = fopen(filename, "r");

        // run parameters
        fgets(buf, 32, fid);
        sscanf(buf, "%lf", &Nphotons); // desired time duration of run [min]
        fgets(buf, 32, fid);
        sscanf(buf, "%lf", &p); // RMT chance of a foward photon doing a bias scattering.
        fgets(buf, 32, fid);
        sscanf(buf, "%d", &Ndetectors); // RMT number of alines
        fgets(buf, 32, fid);
        sscanf(buf, "%lf", &det_radius); // RMT radius of the detector
        fgets(buf, 32, fid);
        sscanf(buf, "%lf", &cos_accept); // RMT cos of the accepted angle where photon is detected
        fgets(buf, 32, fid);
        sscanf(buf, "%d", &Nx); // # of bins
        fgets(buf, 32, fid);
        sscanf(buf, "%d", &Ny); // # of bins
        fgets(buf, 32, fid);
        sscanf(buf, "%d", &Nz); // # of bins

        fgets(buf, 32, fid);
        sscanf(buf, "%f", &dx); // size of bins [cm]
        fgets(buf, 32, fid);
        sscanf(buf, "%f", &dy); // size of bins [cm]
        fgets(buf, 32, fid);
        sscanf(buf, "%f", &dz); // size of bins [cm]
        printf("%lf\n", p);
        // launch parameters
        fgets(buf, 32, fid);
        sscanf(buf, "%f", &radius); // radius

        fgets(buf, 32, fid);
        sscanf(buf, "%f", &zsurf); // z_surface

        // tissue optical properties
        fgets(buf, 32, fid);
        sscanf(buf, "%d", &Nt); // # of tissue types in tissue list

        for (i = 1; i <= Nt; i++)
        {
            fgets(buf, 32, fid);
            sscanf(buf, "%f", &muav[i]); // absorption coeff [cm^-1]
            fgets(buf, 32, fid);
            sscanf(buf, "%f", &musv[i]); // scattering coeff [cm^-1]
            fgets(buf, 32, fid);
            sscanf(buf, "%f", &gv[i]); // anisotropy of scatter [dimensionless]
            fgets(buf, 32, fid);
            sscanf(buf, "%f", &nrv[i]); // refractive index [dimensionless]
            fgets(buf, 32, fid);
            sscanf(buf, "%f", &gf_tt[i]); // forward g for TT
            fgets(buf, 32, fid);
            sscanf(buf, "%f", &gb_tt[i]); // backward g for TT
            fgets(buf, 32, fid);
            sscanf(buf, "%f", &alf_tt[i]); // forward alpha for TT
            fgets(buf, 32, fid);
            sscanf(buf, "%f", &alb_tt[i]); // backward alpha for TT
            fgets(buf, 32, fid);
            sscanf(buf, "%f", &C_tt[i]); // attention
            fgets(buf, 32, fid);
            sscanf(buf, "%f", &eff_r[i]); // attention
            /*TT scattering phase function */
            if(gf_tt[i]!=0){
            for(int j=0;j<nang; j++){
                TT[i][j] = TT_final1(gf_tt[i],gb_tt[i],alf_tt[i],alb_tt[i],C_tt[i],2./(nang-1)*j-1);
                if(TT[i][j]>max_TT_v[i])
                max_TT_v[i]=TT[i][j];
            }
            for(int j=0;j<nang_is; j++){
                TT_is[i][j] = TT_is_final1(gf_tt[i],gb_tt[i],alf_tt[i],alb_tt[i],C_tt[i],1./(nang_is-1)*j);
                if(TT_is[i][j]>max_TT_is_v[i])
                max_TT_is_v[i]=TT_is[i][j];
            }
        }
            /******************/
        }
        fgets(buf, 32, fid);
        sscanf(buf, "%f", &zfocus); // z_focus
        fgets(buf, 32, fid);
        sscanf(buf, "%f", &waist); // waist
        fgets(buf, 32, fid);
        sscanf(buf, "%f", &zr); // zr
        fclose(fid);

        NN = Nx * Ny * Nz;
        Nyx = Nx * Ny;
        v = (char *)malloc(NN * sizeof(char)); /* tissue structure */

        Detsnum = malloc(sizeof(int));

        

        // read binary file
        strcpy(filename, myname);
        strcat(filename, "_T.bin");
        fid = fopen("settings/infi1_T.bin", "rb");
        fread(v, sizeof(char), NN, fid);
        fclose(fid);

        // Show tissue on screen, along central z-axis, by listing tissue type #'s.
        iy = Ny / 2;
        ix = Nx / 2;

        /**************************
         * ============================ MAJOR CYCLE ========================
         **********/

        start_time = clock();

        /**** INITIALIZATIONS
         *****/
        RandomGen(0, -(int)time(NULL) % (1 << 15), NULL);
        Rd = 0.0;

        /**** RUN
         Launch N photons, initializing each one before progation.
         *****/
        printf("------------- Begin Monte Carlo -------------\n");
        printf("%s\n", myname);
        i_photon = 0;
        c_photon = 0;
        c_photon_ss = 0;
        c_photon_ms = 0;
        do
        {
            /**** LAUNCH: Initialize photon position and trajectory *****/
            i_photon += 1; /* increment photon count */
            double dis_rand_old = 0;
            zfocus_used = zfocus;
            back_pick = 0;
            k = 0;
            scatter_find = 0;
            haven_scatter = 0;
            Point3D points_old;
            points_old.x = 1000;
            points_old.y = 1000;
            points_old.z = 1000; 
            int p_old = -1000;
            int p_old_center = -1000;
            double have_skip=0;
            double sleft_max = 0;
            particle_finde_index = -1;
            particle_finde_index_center = -1;

            num_s = 0;
            W = 1.0;               /* set photon weight to one */
            photon_status = ALIVE; /* Launch an ALIVE photon */

            det_num = -1;        /* photon not detected yet */
            first_bias_done = 0; /* photon not biased back - scattered yet */
            cont_exist = 0;      /* no split generated yet */
            L_current = 1;       /* photon 's initial likelihood */
            s_total = 0;         /* photon 's initial path length */
            z_max = 0;           /* photon 's initial depth reached */
            split_num = 0;

            // At 1000th photon, update Nphotons to achieve desired runtime (time_min)
            if (i_photon == 1)
                temp_time = clock();
            if (i_photon == 1000)
            {
                finish_time = clock();
                printf("FUll SIMULATION TIME = %0.2f min for photon numbers = %f \n", Nphotons / Sep_photons * (finish_time - temp_time) / CLOCKS_PER_SEC / 60, Nphotons);
            }

            Pick_det = p_num;
            if (Ndetectors == 1)
            {
                detx = 0;
            }
            else
            {
                detx = 2 * radius * (Pick_det - 1) / (Ndetectors - 1) - radius;
            }

            /**** SET SOURCE****/
            double rnd1 = RandomNum;
            double rnd2 = RandomNum;

            // generateRandomGaussian(&x, &y,rand1[i_photon],rand2[i_photon]);
            generateRandomGaussian(&x, &y,rnd1,rnd2);
            x = detx + det_radius * x;
            y = 0 + det_radius * y;
            z = 0;
            generateGaussianDirection(x, y, z, &ux, &uy, &uz, zr, zfocus_used,detx);
            cos_accept = uz;

            /****************************/
            /* Get tissue voxel properties of launchpoint.
             * If photon beyond outer edge of defined voxels,
             * the tissue equals properties of outermost voxels.
             * Therefore, set outermost voxels to infinite background value.
             */
            ix = floor(Nx / 2 + x / dx);
            iy = floor(Ny / 2 + y / dy);
            iz = floor(z / dz);
            if (ix >= Nx)
                ix = Nx - 1;
            if (iy >= Ny)
                iy = Ny - 1;
            if (iz >= Nz)
                iz = Nz - 1;
            if (ix < 0)
                ix = 0;
            if (iy < 0)
                iy = 0;
            if (iz < 0)
                iz = 0;
            /* Get the tissue type of located voxel */
            i = (long)(iz * Ny * Nx + ix * Ny + iy);
            type = v[i];
            mua = muav[type];
            mus = musv[type];
            g = gv[type];
            nr = nrv[type];
            gf = gf_tt[type];
            gb = gb_tt[type];
            alpha_f = alf_tt[type];
            alpha_b = alb_tt[type];
            C = C_tt[type];
            max_TT = max_TT_v[type];
            max_TT_is = max_TT_is_v[type];
            effr = eff_r[type];

            // initialize as 1 = inside volume, but later check as photon propagates.
            surfflag = 1; // initially inside tissue
            // NOTE: must launch photons at tissue surface, or surfflag will go to 0.
            det_z = zsurf; //

            /* HOP_DROP_SPIN_CHECK
             Propagate one photon until it dies as determined by ROULETTE.
             *******/
            do
            {
                /**** HOP
                 Take step to new position
                 s = dimensionless stepsize
                 x, uy, uz are cosines of current photon trajectory
                 *****/
                // while ((rnd = RandomNum) <= 0.0);              /* yields 0 < rnd <= 1 */
                // sleft = -log(rnd); /* dimensionless step */

                if (photon_status == DEAD)
                { // load the continuing photon and update the flags

                    x = x_cont;
                    y = y_cont;
                    z = z_cont;
                    ux = ux_cont;
                    uy = uy_cont;
                    uz = uz_cont;
                    i = i_cont;
                    s_total = s_total_cont;
                    split_num = split_num_cont;
                    z_max = z_max_cont;
                    num_s = num_s_cont;
                    type = v[i];
                    type_old = type_old_cont;
                    mua = muav[type];
                    mus = musv[type];
                    g = gv[type];
                    nr = nrv[type];
                    gf = gf_tt[type];
                    gb = gb_tt[type];
                    alpha_f = alf_tt[type];
                    alpha_b = alb_tt[type];
                    C = C_tt[type];
                    effr = eff_r[type];
                    max_TT = max_TT_v[type];
                    max_TT_is = max_TT_is_v[type];
                    W = W_cont;
                    L_current = L_cont;
                    cont_exist = 0;
                    photon_status = ALIVE;
                    first_bias_done = 0;
                    det_num = -1;
                    have_skip=1;
                }

                /*particle fixed calculation*/
                if(type!=4 & type!=5){
                    while ((rnd = RandomNum) <= 0.0);              /* yields 0 < rnd <= 1 */
                    sleft = -log(rnd); /* dimensionless step */
                    if(type!=1)
                    scatter_find=1;
                }
                else if(type==4){
                Line3D line;
                line.start.x = x;line.start.y = y ;line.start.z = z ;
                ix_particle = floor(Nx / 2 + line.start.x / dx);
                iy_particle = floor(Ny / 2 + line.start.y / dy);
                iz_particle = floor(line.start.z / dz);
                if(ix_particle<0){
                    ix_particle = 0;
                    line.start.x = -Nx / 2*dx;
                }
                if(iy_particle<0){
                    iy_particle = 0;
                    line.start.y = -Ny / 2*dy;
                }
                if(iz_particle<0){
                    iz_particle = 0;
                    line.start.z = 0;
                }
                if(ix_particle>=Nx){
                    ix_particle = Nx-1;
                    line.start.x = ((Nx-1)-(Nx / 2))*dx;
                }
                if(iy_particle>=Ny){
                    iy_particle = Ny-1;
                    line.start.y = ((Ny-1)-(Ny / 2))*dy;
                }
                if(iz_particle>=Nz){
                    iz_particle = Nz-1;
                    line.start.z = (Nz-1)*dz;
                }
                i_particle = (long)(iz_particle * Ny * Nx + ix_particle * Ny + iy_particle);
                double minX,maxX,minY,maxY,minZ,maxZ;
                int lowerIndex_cur = -1;
                int upperIndex_cur = -1;
                int lowerIndex_next = -1;
                int upperIndex_next = -1;
                findIndices(particle_index, particle_num, i_particle-1, &lowerIndex_cur, &upperIndex_cur);
                sleft_max = ls+FindVoxelFace2(&fdir, line.start.x, line.start.y, line.start.z, &det_num, Pick_det, detx, det_radius, det_z, cos_accept, Ndetectors, dx, dy, dz, ux, uy, uz);             
                if(upperIndex_cur<0){
                    sleft = sleft_max;
                }
                else{
                    findIndices(particle_index, particle_num, i_particle, &lowerIndex_next, &upperIndex_next);
                    if(lowerIndex_next<0)
                    sleft = sleft_max;
                    else{
                    line.end.x = x+sleft_max*ux;
                    line.end.y = y+sleft_max*uy;
                    line.end.z = z+sleft_max*uz;
                    scatter_dis = 1000;
                    sleft = sleft_max;
                    for(int p1 = lowerIndex_cur+1;p1<=lowerIndex_next;p1++){
                        Point3D points;
                        points.x = particles[p1][0];
                        points.y = particles[p1][1];
                        points.z = particles[p1][2];
                        double distance = distanceToLineSegment(points, line.start, line.end);
                        if(distance<=effr){
                            double dis_temp = sqrt((points.x-line.start.x)*(points.x-line.start.x)+(points.y-line.start.y)*(points.y-line.start.y)+(points.z-line.start.z)*(points.z-line.start.z));
                            double dis_sleft = sqrt(dis_temp*dis_temp-distance*distance);
                            double dis_rand = sqrt(effr*effr-distance*distance);
                            dis_sleft = dis_sleft+dis_rand+1e-9;
                            double x_temp = x+ux*dis_sleft;
                            double y_temp = y+uy*dis_sleft;
                            double z_temp = z+uz*dis_sleft;
                            double dis_temp2 = sqrt((points.x-x_temp)*(points.x-x_temp)+(points.y-y_temp)*(points.y-y_temp)+(points.z-z_temp)*(points.z-z_temp));
                            if(dis_temp<scatter_dis&fabs(dis_temp2-effr)<1e-5){
                                scatter_dis = dis_temp;
                                sleft = dis_sleft;
                                scatter_find=1;
                                particle_finde_index = p1;
                            } 
                        }
                    }
                    if(have_skip==0&particle_finde_index>=0&first_bias_done == 1&uz<0){
                        if(p_old==particle_finde_index){
                            num_s-=1;
                            particle_finde_index=-1;
                            have_skip=1;
                            p_old=-1000;  
                        }
                    }
                    if(particle_finde_index>=0){
                        p_old = particle_finde_index;
                    }
                }
                }
                }
                else{
                    Line3D line;
                line.start.x = x;line.start.y = y ;line.start.z = z ;
                ix_particle = floor(Nx / 2 + line.start.x / dx);
                iy_particle = floor(Ny / 2 + line.start.y / dy);
                iz_particle = floor(line.start.z / dz);
                if(ix_particle<0){
                    ix_particle = 0;
                    line.start.x = -Nx / 2*dx;
                }
                if(iy_particle<0){
                    iy_particle = 0;
                    line.start.y = -Ny / 2*dy;
                }
                if(iz_particle<0){
                    iz_particle = 0;
                    line.start.z = 0;
                }
                if(ix_particle>=Nx){
                    ix_particle = Nx-1;
                    line.start.x = ((Nx-1)-(Nx / 2))*dx;
                }
                if(iy_particle>=Ny){
                    iy_particle = Ny-1;
                    line.start.y = ((Ny-1)-(Ny / 2))*dy;
                }
                if(iz_particle>=Nz){
                    iz_particle = Nz-1;
                    line.start.z = (Nz-1)*dz;
                }
                i_particle_center = (long)(iz_particle * Ny * Nx + ix_particle * Ny + iy_particle);
                double minX,maxX,minY,maxY,minZ,maxZ;
                int lowerIndex_cur = -1;
                int upperIndex_cur = -1;
                int lowerIndex_next = -1;
                int upperIndex_next = -1;
                findIndices(particle_center_index, particle_center_num, i_particle_center-1, &lowerIndex_cur, &upperIndex_cur);
                sleft_max = ls+FindVoxelFace2(&fdir, line.start.x, line.start.y, line.start.z, &det_num, Pick_det, detx, det_radius, det_z, cos_accept, Ndetectors, dx, dy, dz, ux, uy, uz);              
                if(upperIndex_cur<0){
                    sleft = sleft_max;
                }
                else{
                    findIndices(particle_center_index, particle_center_num, i_particle_center, &lowerIndex_next, &upperIndex_next);
                    if(lowerIndex_next<0)
                    sleft = sleft_max;
                    else{
                    line.end.x = x+sleft_max*ux;
                    line.end.y = y+sleft_max*uy;
                    line.end.z = z+sleft_max*uz;
                    scatter_dis = 1000;
                    sleft = sleft_max;
                    for(int p1 = lowerIndex_cur+1;p1<=lowerIndex_next;p1++){
                        Point3D points;
                        points.x = particles_center[p1][0];
                        points.y = particles_center[p1][1];
                        points.z = particles_center[p1][2];
                        double distance = distanceToLineSegment(points, line.start, line.end);
                        if(distance<=effr){
                            double dis_temp = sqrt((points.x-line.start.x)*(points.x-line.start.x)+(points.y-line.start.y)*(points.y-line.start.y)+(points.z-line.start.z)*(points.z-line.start.z));
                            double dis_sleft = sqrt(dis_temp*dis_temp-distance*distance);
                            double dis_rand = sqrt(effr*effr-distance*distance);
                            dis_sleft = dis_sleft+dis_rand+1e-9;
                            double x_temp = x+ux*dis_sleft;
                            double y_temp = y+uy*dis_sleft;
                            double z_temp = z+uz*dis_sleft;
                            double dis_temp2 = sqrt((points.x-x_temp)*(points.x-x_temp)+(points.y-y_temp)*(points.y-y_temp)+(points.z-z_temp)*(points.z-z_temp));
                            if(dis_temp<scatter_dis&fabs(dis_temp2-effr)<1e-5){
                                scatter_dis = dis_temp;
                                sleft = dis_sleft;
                                scatter_find=1;
                                particle_finde_index_center = p1;
                            } 
                        }
                    }
                    if(have_skip==0&particle_finde_index_center>=0&first_bias_done == 1&uz<0){
                        if(p_old_center==particle_finde_index_center){
                            num_s-=1;
                            particle_finde_index_center=-1;
                            have_skip=1;
                            p_old_center=-1000; 
                        }
                    }
                    if(particle_finde_index_center>=0){
                        p_old_center = particle_finde_index_center;
                    }
                }
                }
                }
            
                do
                {                            // while sleft>0
                    if(type!=4 & type!=5)
                        s = sleft / (mus); /* Step size [cm].*/
                    else{                         // while sleft>0
                    s = sleft / 1; /* Step size [cm].*/
                    }

                    if(haven_scatter==0){
                        double uz_now = uz;
                        generateGaussianDirection(x, y, z, &ux, &uy, &uz, zr, zfocus_used,detx);
                        if(uz_now>0) uz=uz;
                        else uz=-uz;
                    }
                    tempx = x + s * ux; /* Update positions. [cm] */
                    tempy = y + s * uy;
                    tempz = z + s * uz;

                    sv = SameVoxel(x, y, z, tempx, tempy, tempz, dx, dy, dz);
                    if (sv) /* photon in same voxel */
                    {
                        x = tempx; /* Update positions. */
                        y = tempy;
                        z = tempz;
                        s_total += s*nr;
                        // }
                        /**** DROP
                         Drop photon weight (W) into local bin.
                         *****/
                        if(type != 4 & type != 5){
                            absorb = W * (1 - exp(-mua * s));
                        /* photon weight absorbed at this step */
                        }
                        else{
                           absorb = (mua/(mua+mus)) *W;
                        }
                        W -= absorb;
                        /* Update sleft */
                        sleft = 0; /* dimensionless step remaining */
                        type_old = type;
                    }
                    else /* photon has crossed voxel boundary */
                    {
                        /* step to voxel face + "littlest step" so just inside new voxel. */
                        s = ls + FindVoxelFace2(&fdir, x, y, z, &det_num, Pick_det, detx, det_radius, det_z, cos_accept, Ndetectors, dx, dy, dz, ux, uy, uz);

                        /*** DROP: Drop photon weight (W) into local bin  ***/
                        if(type != 4 & type != 5){
                            absorb = W * (1 - exp(-mua * s));
                        /* photon weight absorbed at this step */
                        }
                        else{
                           absorb = mua/(mua+mus) *W;
                        }
                        W -= absorb;
                    
                        if (det_num != -1)
                        { /* check if the photon is detected . */

                            /* Update total path length */
                            s_total += s*nr;
                            x = x + s * ux; /* Update positions. [cm] */
                            y = y + s * uy;
                            z = z + s * uz;
                            s_total_sc = sqrt((x-x_last)*(x-x_last)+(y-y_last)*(y-y_last)+(z-z_last)*(z-z_last));
                          
							
                            /* Save properties of interest */
                            if (L_current > 0 & L_current <= 1 & det_num == Pick_det)
                            { // avoid NAN and zero likelihood, and avoid cross - detection
                                if (aa == 1)
                                {
                                    DetS = realloc(DetS, (c_photon + 1) * sizeof(float));
                                    DetS[c_photon] = s_total;
                                    DetL = realloc(DetL, (c_photon + 1) * sizeof(float));
                                    DetL[c_photon] = L_current;
                                    DetW = realloc(DetW, (c_photon + 1) * sizeof(float));
                                    DetW[c_photon] = W;

                                }
                                else
                                {
                                    DetS = realloc(DetS, (c_photon_mask[aa - 2] + c_photon + 1) * sizeof(float));
                                    DetS[c_photon_mask[aa - 2] + c_photon] = s_total;
                                    DetL = realloc(DetL, (c_photon_mask[aa - 2] + c_photon + 1) * sizeof(float));
                                    DetL[c_photon_mask[aa - 2] + c_photon] = L_current;
                                    DetW = realloc(DetW, (c_photon_mask[aa - 2] + c_photon + 1) * sizeof(float));
                                    DetW[c_photon_mask[aa - 2] + c_photon] = W;
                                }
                                
                                /* increment collected photon count */
                                Detsnum = realloc(Detsnum, (c_photon + 1) * sizeof(int));
                                Detsnum[c_photon] = num_s;
                                c_photon += 1;
                            }

                            if (L_current > 0 & L_current <= 1 & det_num == Pick_det & num_s == 1)
                            { // avoid NAN and zero likelihood, and avoid cross - detection
                                if (aa == 1)
                                {
                                    DetS_ss = realloc(DetS_ss, (c_photon_ss + 1) * sizeof(float));
                                    DetS_ss[c_photon_ss] = s_total;
                                    DetL_ss = realloc(DetL_ss, (c_photon_ss + 1) * sizeof(float));
                                    DetL_ss[c_photon_ss] = L_current;
                                    DetW_ss = realloc(DetW_ss, (c_photon_ss + 1) * sizeof(float));
                                    DetW_ss[c_photon_ss] = W;

                                }
                                else
                                {
                                    DetS_ss = realloc(DetS_ss, (c_photon_mask_ss[aa - 2] + c_photon_ss + 1) * sizeof(float));
                                    DetS_ss[c_photon_mask_ss[aa - 2] + c_photon_ss] = s_total;
                                    DetL_ss = realloc(DetL_ss, (c_photon_mask_ss[aa - 2] + c_photon_ss + 1) * sizeof(float));
                                    DetL_ss[c_photon_mask_ss[aa - 2] + c_photon_ss] = L_current;
                                    DetW_ss = realloc(DetW_ss, (c_photon_mask_ss[aa - 2] + c_photon_ss + 1) * sizeof(float));
                                    DetW_ss[c_photon_mask_ss[aa - 2] + c_photon_ss] = W;
                                }
                                
                                /* increment collected photon count */
                                c_photon_ss += 1;
                            }

                            if (L_current > 0 & L_current <= 1 & det_num == Pick_det& num_s > 1)
                            { // avoid NAN and zero likelihood, and avoid cross - detection
                                if (aa == 1)
                                {
                                    DetS_ms = realloc(DetS_ms, (c_photon_ms + 1) * sizeof(float));
                                    DetS_ms[c_photon_ms] = s_total;
                                    DetL_ms = realloc(DetL_ms, (c_photon_ms + 1) * sizeof(float));
                                    DetL_ms[c_photon_ms] = L_current;
                                    DetW_ms = realloc(DetW_ms, (c_photon_ms + 1) * sizeof(float));
                                    DetW_ms[c_photon_ms] = W;

                                }
                                else
                                {
                                    DetS_ms = realloc(DetS_ms, (c_photon_mask_ms[aa - 2] + c_photon_ms + 1) * sizeof(float));
                                    DetS_ms[c_photon_mask_ms[aa - 2] + c_photon_ms] = s_total;
                                    DetL_ms = realloc(DetL_ms, (c_photon_mask_ms[aa - 2] + c_photon_ms + 1) * sizeof(float));
                                    DetL_ms[c_photon_mask_ms[aa - 2] + c_photon_ms] = L_current;
                                    DetW_ms = realloc(DetW_ms, (c_photon_mask_ms[aa - 2] + c_photon_ms + 1) * sizeof(float));
                                    DetW_ms[c_photon_mask_ms[aa - 2] + c_photon_ms] = W;
                                }
                                
                                /* increment collected photon count */
                                c_photon_ms += 1;
                            }
                            photon_status = DEAD;
                            sleft = 0;
                        }
                        else
                        {
                            /* Update sleft */
                            if(type!=4 & type != 5){
                                sleft -= s * (mus);
                            }
                            else{
                                sleft -= s * 1; /* dimensionless step remaining */
                            }
                            if (sleft <= ls)
                                sleft = 0;

                            /* Update positions. */
                            x += s * ux;
                            y += s * uy;
                            z += s * uz;

                            /* Update total path length */ // 
                            s_total += s*nr;
                            

                            // pointers to voxel containing optical properties
                            ix = floor(Nx / 2 + x / dx);
                            iy = floor(Ny / 2 + y / dy);
                            iz = floor(z / dz);

                            //*** ESCAPE or not
                            if ((surfflag == 1) & (z <= zsurf)) // escape at surface
                            {
                                Rd += W;
                                i = (long)(Nx * ix + iy);
                                surfflag = 0; // disable repeated assignments to Rd, R[i]
                            }
                            if (z < 0) // escape cube
                            {
                                photon_status = DEAD;
                                sleft = 0;
                            }
                            else // No escape
                            {

                                if (iz >= Nz)
                                {
                                    iz = Nz - 1;
                                    photon_status = DEAD;
                                    sleft = 0;
                                }
                                if (ix >= Nx)
                                {
                                    ix = Nx - 1;
                                    photon_status = DEAD;
                                    sleft = 0;
                                }
                                if (iy >= Ny)
                                {
                                    iy = Ny - 1;
                                    photon_status = DEAD;
                                    sleft = 0;
                                }
                                if (iz < 0)
                                {
                                    iz = 0;
                                    photon_status = DEAD;
                                    sleft = 0;
                                }
                                if (ix < 0)
                                {
                                    ix = 0;
                                    photon_status = DEAD;
                                    sleft = 0;
                                }
                                if (iy < 0)
                                {
                                    iy = 0;
                                    photon_status = DEAD;
                                    sleft = 0;
                                }

                                type_old = type;
                                double nr_old = nrv[type_old];
                                // update pointer to tissue type
                                i = (long)(iz * Ny * Nx + ix * Ny + iy);
                                type = v[i];
                                mua = muav[type];
                                mus = musv[type];
                                g = gv[type];
                                nr = nrv[type];
                                gf = gf_tt[type];
                                gb = gb_tt[type];
                                alpha_f = alf_tt[type];
                                alpha_b = alb_tt[type];
                                C = C_tt[type];
                                max_TT = max_TT_v[type];
                                max_TT_is = max_TT_is_v[type];
                                effr = eff_r[type];

                            /*Reflection and refraction, NOT accurate for the voxel-based model, so remove it*/
                            //     if(photon_status == ALIVE){
                            //     if(fdir == 3){
                            //         if(uz < 0.0){
                            //         if(nr != nr_old)
                            //             CrossUpOrNot(nr_old,nr);
                            //     }
                                    
                            //     else{
                            //         if(nr != nr_old)
                            //             CrossDnOrNot(nr_old,nr);
                            //     }
                            //     }
                            //     else if(fdir == 2){
                            //         if(uy < 0.0){
                            //         if(nr != nr_old)
                            //             CrossOutOrNot(nr_old,nr);
                            //     }
                                    
                            //     else{
                            //         if(nr != nr_old)
                            //             CrossInOrNot(nr_old,nr);
                            //     }
                            //     }
                            //     else{
                            //         if(ux < 0.0){
                            //         if(nr != nr_old)
                            //             CrossLeftOrNot(nr_old,nr);
                            //     }
                                    
                            //     else{
                            //         if(nr != nr_old)
                            //             CrossRightOrNot(nr_old,nr);
                            //     }
                            //     }
                            // }
                            // i = (long)(iz * Ny * Nx + ix * Ny + iy);
                            // type = v[i];
                            // mua = muav[type];
                            // mus = musv[type];
                            // g = gv[type];
                            // nr = nrv[type];
                            // gf = gf_tt[type];
                            // gb = gb_tt[type];
                            // alpha_f = alf_tt[type];
                            // alpha_b = alb_tt[type];
                            // C = C_tt[type];
                            // max_TT = max_TT_v[type];
                            // max_TT_is = max_TT_is_v[type];
                            // effr = eff_r[type];

                                if(type_old!=4 & type_old!=5){
                                    if(type==4 || type==5){
                                        sleft=0;
                                        sleft_max=sleft;
                                        scatter_find=0;
                                        if(haven_scatter==0){
                                            double uz_now = uz;
                                            generateGaussianDirection(x, y, z, &ux, &uy, &uz, zr, zfocus_used,detx);
                                            if(uz_now>0) uz=uz;
                                            else uz=-uz;
                                        }
                                    }
                                }
                            }
                        }
                    } //(sv) /* same voxel */

                } while (sleft > 0); // do...while
                /**** SPIN AND SPLIT*/

                if (photon_status == ALIVE && g != 1&&scatter_find==1)
                {
                    x_last = x;
                    y_last = y;
                    z_last = z;
                    if (first_bias_done == 0 & uz > 0)
                    { /* apply the first biased scattering */
                        /* Sample for costheta_B */
                        rnd = RandomNum;
                        find_sample = 0;
                        if (g == 0.0)
                            costheta_B = 2.0 * rnd - 1.0;
                        else
                        {
                            if (g != 1)
                            {
                                costheta_B = cos_is_final1(gf,gb,alpha_f,alpha_b,C);
                            }
                            else
                            {
                                costheta_B = 1;
                            }
                        }

                        sintheta_B = sqrt(1.0 - costheta_B * costheta_B);
                        /* Sample psi . */
                        psi = 2.0 * PI * RandomNum;
                        cospsi = cos(psi);
                        if (psi < PI)
                            sinpsi = sqrt(1.0 - cospsi * cospsi); /* sqrt () is faster than sin (). */
                        else
                            sinpsi = -sqrt(1.0 - cospsi * cospsi);
                        /* Compute the unit vector v towards the actual position of the detector , ...
                         where detx is chosen uniformly along the centers of the collecting fiber ...
                         array . */
                        if (Ndetectors == 1)
                        {
                            detx = 0;
                        }
                        else
                        {
                            detx = 2 * radius * (Pick_det - 1) / (Ndetectors - 1) - radius;
                        }

                        dety = 0;
                        detz = det_z;
                        temp = sqrt((x - detx) * (x - detx) + (y - dety) * (y - dety) + (z - detz) * (z - detz));
                        vx = -(x - detx) / temp;
                        vy = -(y - dety) / temp;
                        vz = -(z - detz) / temp;

                        /* New trajectory u' = (upx , upy , upz) */
                        if (1 - fabs(vz) <= ONE_MINUS_COSZERO)
                        { /* close to perpendicular . */
                            upx = sintheta_B * cospsi;
                            upy = sintheta_B * sinpsi;
                            upz = costheta_B * SIGN(vz); /* SIGN () is faster than division . */
                        }
                        else
                        { /* usually use this option */
                            temp = sqrt(1.0 - vz * vz);
                            upx = sintheta_B * (vx * vz * cospsi - vy * sinpsi) / temp + vx * costheta_B;
                            upy = sintheta_B * (vy * vz * cospsi + vx * sinpsi) / temp + vy * costheta_B;
                            upz = -sintheta_B * cospsi * temp + vz * costheta_B;
                        }
                        /* Compute the likelihood ratio for this particular biased ...
                         back - scattering */

                        costheta_S = upx * ux + upy * uy + upz * uz;
                        
                        temp11 = TT_final1(gf, gb, alpha_f, alpha_b, C, costheta_S);
                        temp22 = TT_is_final1(gf, gb, alpha_f, alpha_b, C, costheta_B);
                        L_temp = temp11 / temp22;

                        if (L_temp < 1 * (1 - ls) && split_num > -1)
                        { // yes , do the unbiased spin and save the trajectory for the continuing photon packet
                            L_cont = L_current * (1 - L_temp);
                            i_cont = i;
                            split_num += 1;
                            /* the unbiased spin */
                            /* Sample for costheta */
                            rnd = RandomNum;
                            if (g == 0.0)
                                costheta = 2.0 * rnd - 1.0;
                            else
                            {
                                find_sample = 0;

                                if (g != 1)
                                {
                                    while (costheta < 0)
                                    {
                                        costheta = cos_final1(gf,gb,alpha_f,alpha_b,C);
                                    }
                                }
                                else
                                {
                                    costheta = 1;
                                }
                            }
                            sintheta = sqrt(1.0 - costheta * costheta); /* sqrt () is faster than sin (). */
                            /* Sample psi . */

                            psi = 2.0 * PI * RandomNum;
                            cospsi = cos(psi);
                            if (psi < PI)
                                sinpsi = sqrt(1.0 - cospsi * cospsi); /* sqrt () is faster than sin (). */
                            else
                                sinpsi = -sqrt(1.0 - cospsi * cospsi);
                            /* New trajectory . */

                            if (1 - fabs(uz) <= ONE_MINUS_COSZERO)
                            { /* close to perpendicular . */
                                uxx = sintheta * cospsi;
                                uyy = sintheta * sinpsi;
                                uzz = costheta * SIGN(uz); /* SIGN () is faster than division . */
                            }
                            else
                            { /* usually use this option */

                                temp = sqrt(1.0 - uz * uz);
                                uxx = sintheta * (ux * uz * cospsi - uy * sinpsi) / temp + ux * costheta;
                                uyy = sintheta * (uy * uz * cospsi + ux * sinpsi) / temp + uy * costheta;
                                uzz = -sintheta * cospsi * temp + uz * costheta;
                            }
                            ux_cont = uxx;
                            uy_cont = uyy;
                            uz_cont = uzz;

                            x_cont = x;
                            y_cont = y;
                            z_cont = z;
                            W_cont = W;
                            s_total_cont = s_total;
                            z_max_cont = z_max;
                            num_s += 1;
                            num_s_cont = num_s;
                            split_num_cont = split_num;
                            L_current *= L_temp;
                            cont_exist = 1;
                            type_old_cont = type_old;
                        }
                        else
                        { // no continuing photon packet
                            L_current *= L_temp;
                            cont_exist = 0;
                            num_s += 1;
                        }
                        /* Update trajectory */
                        ux = upx;
                        uy = upy;
                        uz = upz;
                        first_bias_done = 1;
                        scatter_find=0;
                        haven_scatter = 1;
                        particle_finde_index = -1;
                        particle_finde_index_center = -1;
                    }
                    else
                    { /* first biased back - scattering already done , apply additional biased ...
                   forward - scattering */
                        if (RandomNum <= p )
                        { // apply biased forward - scattering
                            /* Sample for costheta_B */
                            rnd = RandomNum;
                            find_sample = 0;
                            if (g == 0.0)
                                costheta_B = 2.0 * rnd - 1.0;
                            else
                            {

                                if (g != 1)
                                {
                                    costheta_B = cos_is_final1(gf,gb,alpha_f,alpha_b,C);
                                }
                                else
                                {
                                    costheta_B = 1;
                                }
                            }
                            sintheta_B = sqrt(1.0 - costheta_B * costheta_B);
                            /* Sample psi . */
                            psi = 2.0 * PI * RandomNum;
                            cospsi = cos(psi);
                            if (psi < PI)
                                sinpsi = sqrt(1.0 - cospsi * cospsi); /* sqrt () is faster than sin (). */
                            else
                                sinpsi = -sqrt(1.0 - cospsi * cospsi);
                            /* Compute the unit vector v towards the actual position of the ...
                            detector , where detx is chosen uniformly along the centers of the ...
                            collecting fiber array . */
                            if (Ndetectors == 1)
                                detx = 0;
                            else
                                detx = 2 * radius * (Pick_det - 1) / (Ndetectors - 1) - radius;
                            dety = 0;
                            detz = det_z;
                            temp = sqrt((x - detx) * (x - detx) + (y - dety) * (y - dety) + (z - detz) * (z - detz));
                            vx = -(x - detx) / temp;
                            vy = -(y - dety) / temp;
                            vz = -(z - detz) / temp;
                            /* New trajectory u' = (upx , upy , upz) */
                            if (1 - fabs(vz) <= ONE_MINUS_COSZERO)
                            { /* close to perpendicular . */
                                upx = sintheta_B * cospsi;
                                upy = sintheta_B * sinpsi;
                                upz = costheta_B * SIGN(vz); /* SIGN () is faster than division . */
                            }
                            else
                            { /* usually use this option */
                                temp = sqrt(1.0 - vz * vz);
                                upx = sintheta_B * (vx * vz * cospsi - vy * sinpsi) / temp + vx * costheta_B;
                                upy = sintheta_B * (vy * vz * cospsi + vx * sinpsi) / temp + vy * costheta_B;
                                upz = -sintheta_B * cospsi * temp + vz * costheta_B;
                            }
                            /* Compute the likelihood ratio for this particular biased ...
                           forward - scattering */
                            costheta_S = upx * ux + upy * uy + upz * uz;
                            
                            f_TT = TT_final1(gf, gb, alpha_f, alpha_b, C, costheta_S);
                            f_B = TT_is_final1(gf, gb, alpha_f, alpha_b, C, costheta_B);
                            L_temp = f_TT / (p * f_B + (1 - p) * f_TT);
                            L_current *= L_temp;
                            /* Update trajectory */
                            ux = upx;
                            uy = upy;
                            uz = upz;
                            num_s += 1;
                            haven_scatter = 1;
                        }
                        else
                        { // apply unbiased scattering
                            /* Sample for costheta */
                            rnd = RandomNum;
                            if (g == 0.0)
                                costheta = 2.0 * rnd - 1.0;
                            else
                            {
                                find_sample = 0;

                                if (g != 1)
                                {
                                    costheta = cos_final1(gf,gb,alpha_f,alpha_b,C);
                                }
                                else
                                {
                                    costheta = 1;
                                }
                            }
                            sintheta = sqrt(1.0 - costheta * costheta); /* sqrt () is faster than sin (). */
                            /* Sample psi . */
                            psi = 2.0 * PI * RandomNum;
                            cospsi = cos(psi);
                            if (psi < PI)
                                sinpsi = sqrt(1.0 - cospsi * cospsi); /* sqrt () is faster than sin (). */
                            else
                                sinpsi = -sqrt(1.0 - cospsi * cospsi);
                            /* New trajectory . */
                            if (1 - fabs(uz) <= ONE_MINUS_COSZERO)
                            { /* close to perpendicular . */
                                uxx = sintheta * cospsi;
                                uyy = sintheta * sinpsi;
                                uzz = costheta * SIGN(uz); /* SIGN () is faster than division . */
                            }
                            else
                            { /* usually use this option */
                                temp = sqrt(1.0 - uz * uz);
                                uxx = sintheta * (ux * uz * cospsi - uy * sinpsi) / temp + ux * costheta;
                                uyy = sintheta * (uy * uz * cospsi + ux * sinpsi) / temp + uy * costheta;
                                uzz = -sintheta * cospsi * temp + uz * costheta;
                            }
                            /* Compute the unit vector v towards the actual position of the ...
                            detector , where detx is chosen uniformly along the centers of the ...
                            collecting fiber array . */
                            if (Ndetectors == 1)
                                detx = 0;
                            else
                                detx = 2 * radius * (Pick_det - 1) / (Ndetectors - 1) - radius;
                            dety = 0;
                            detz = det_z;
                            temp = sqrt((x - detx) * (x - detx) + (y - dety) * (y - dety) + (z - detz) * (z - detz));
                            vx = -(x - detx) / temp;
                            vy = -(y - dety) / temp;
                            vz = -(z - detz) / temp;
                            /* Compute the likelihood ratio for this particular unbiased ...
                            forward - scattering */
                            costheta_S = costheta;
                            costheta_B = uxx * vx + uyy * vy + uzz * vz;
                            
                            f_TT = TT_final1(gf, gb, alpha_f, alpha_b, C, costheta_S);
                            f_B = TT_is_final1(gf, gb, alpha_f, alpha_b, C, costheta_B);
                            L_temp = f_TT / (p * f_B + (1 - p) * f_TT);
                            L_current *= L_temp;
                            /* Update trajectory */
                            ux = uxx;
                            uy = uyy;
                            uz = uzz;
                            num_s += 1;
                            haven_scatter = 1;
                        }
                        scatter_find=0;
                        particle_finde_index = -1;
                        particle_finde_index_center = -1;
                    }

                    /**** CHECK ROULETTE    
               *****/

                    if (W < THRESHOLD)
                    {
                        if (RandomNum <= CHANCE)
                            W /= CHANCE;
                        else
                            photon_status = DEAD;
                    }
                }
            } while (photon_status == ALIVE || cont_exist == 1); /* end STEP_CHECK_HOP_SPIN */

        } while (i_photon < Nphotons); /* end RUN */
        printf("collected photons = %ld\n", c_photon);
        printf("collected SS photons = %ld\n", c_photon_ss);
        printf("collected MS photons = %ld\n", c_photon_ms);
        if (aa == 1){
            c_photon_mask[aa - 1] = c_photon;
            c_photon_mask_ss[aa - 1] = c_photon_ss;
            c_photon_mask_ms[aa - 1] = c_photon_ms;    
        }
            
        else{
            c_photon_mask[aa - 1] = c_photon_mask[aa - 2] + c_photon;
            c_photon_mask_ss[aa - 1] = c_photon_mask_ss[aa - 2] + c_photon_ss;
            c_photon_mask_ms[aa - 1] = c_photon_mask_ms[aa - 2] + c_photon_ms;
        }
            

        printf("------------------------------------------------------\n");
        finish_time = clock();
        time_min = (double)(finish_time - start_time) / CLOCKS_PER_SEC / 60;
        printf("Elapsed Time for %0.3e photons = %5.3f min\n", Nphotons, time_min);
        printf("%0.2e photons per minute\n", Nphotons / time_min);

        /**** process
         Process data
         *****/

        linspace(lambda_start, lambda_end, samplePoints, lambda);
        linspace(2 * PI / lambda_start, 2 * PI / lambda_end, samplePoints, k_lin);
        linspace(1, samplePoints, samplePoints, temp_x);
        for (int i = 0; i < samplePoints; i++)
        {
            k_nolin[i] = 2 * PI / lambda[i];
        }
        z_min = 2 * PI / (k_lin[1] - k_lin[2]) / samplePoints;
        if (aa == 1)
        {
            for (temp_lin = 0; temp_lin < c_photon; temp_lin++)
            {
                for (temp_lin2 = 0; temp_lin2 < samplePoints; temp_lin2++)
                {
                    if (DetS[temp_lin] >= (0 + (temp_lin2 * z_min)) && DetS[temp_lin] < (z_min + (temp_lin2 * z_min)))
                    {
                        MM[temp_lin2][aa - 1] += DetW[temp_lin]*DetL[temp_lin];
                    }
                }
            }
            for (temp_lin = 0; temp_lin < c_photon_ss; temp_lin++)
            {
                for (temp_lin2 = 0; temp_lin2 < samplePoints; temp_lin2++)
                {
                    if (DetS_ss[temp_lin] >= (0 + (temp_lin2 * z_min)) && DetS_ss[temp_lin] < (z_min + (temp_lin2 * z_min)))
                    {
                        MM_ss[temp_lin2][aa - 1] += DetW_ss[temp_lin]*DetL_ss[temp_lin];
                    }
                }
            }
            for (temp_lin = 0; temp_lin < c_photon_ms; temp_lin++)
            {
                for (temp_lin2 = 0; temp_lin2 < samplePoints; temp_lin2++)
                {
                    if (DetS_ms[temp_lin] >= (0 + (temp_lin2 * z_min)) && DetS_ms[temp_lin] < (z_min + (temp_lin2 * z_min)))
                    {
                        MM_ms[temp_lin2][aa - 1] += DetW_ms[temp_lin]*DetL_ms[temp_lin];
                    }
                }
            }
        }
        else
        {
            for (temp_lin = c_photon_mask[aa - 2]; temp_lin < c_photon_mask[aa - 2] + c_photon; temp_lin++)
            {
                for (temp_lin2 = 0; temp_lin2 < samplePoints; temp_lin2++)
                {
                    if (DetS[temp_lin] >= (0 + (temp_lin2 * z_min)) && DetS[temp_lin] < (z_min + (temp_lin2 * z_min)))
                    {
                        MM[temp_lin2][aa - 1] += DetW[temp_lin]*DetL[temp_lin];
                    }
                }
            }
            for (temp_lin = c_photon_mask_ss[aa - 2]; temp_lin < c_photon_mask_ss[aa - 2] + c_photon_ss; temp_lin++)
            {
                for (temp_lin2 = 0; temp_lin2 < samplePoints; temp_lin2++)
                {
                    if (DetS_ss[temp_lin] >= (0 + (temp_lin2 * z_min)) && DetS_ss[temp_lin] < (z_min + (temp_lin2 * z_min)))
                    {
                        MM_ss[temp_lin2][aa - 1] += DetW_ss[temp_lin]*DetL_ss[temp_lin];
                    }
                }
            }
            for (temp_lin = c_photon_mask_ms[aa - 2]; temp_lin < c_photon_mask_ms[aa - 2] + c_photon_ms; temp_lin++)
            {
                for (temp_lin2 = 0; temp_lin2 < samplePoints; temp_lin2++)
                {
                    if (DetS_ms[temp_lin] >= (0 + (temp_lin2 * z_min)) && DetS_ms[temp_lin] < (z_min + (temp_lin2 * z_min)))
                    {
                        MM_ms[temp_lin2][aa - 1] += DetW_ms[temp_lin]*DetL_ms[temp_lin];
                    }
                }
            }
        }

        /**** process
         Process data
         *****/
        printf("%s is done.\n", myname);

        printf("------------------------------------------------------\n");
        now = time(NULL);
        printf("%s\n", ctime(&now));

        free(Detsnum);
        free(v);
        
    }
    end_time_all = clock();
    time_min = (double)(end_time_all - start_time_all) / CLOCKS_PER_SEC / 60;
    printf("Elapsed Time  = %5.3f min\n", time_min);

    // start_time_post = clock();

    /*all photons*/
    linspace(1, samplePoints, samplePoints, temp_x);

    for (temp_lin = 0; temp_lin < samplePoints; temp_lin++)
    {
        sum_fit = 0;
        sum_MM = 0;
        fit_polynomial(temp_x, MM[temp_lin], samplePoints, coeffs);
        for (temp_lin2 = 0; temp_lin2 < samplePoints; temp_lin2++)
        {
            fit_temp[temp_lin2] = polyval(coeffs, ORDER, temp_x[temp_lin2]);
            if (fit_temp[temp_lin2] < 0)
                fit_temp[temp_lin2] = 0;
            sum_fit += fit_temp[temp_lin2];
            sum_MM += MM[temp_lin][temp_lin2];
        }
        for (temp_lin2 = 0; temp_lin2 < samplePoints; temp_lin2++)
        {
            fit_temp[temp_lin2] = fit_temp[temp_lin2] / (sum_fit + ls);
            fit_temp_MM[temp_lin2] = fit_temp[temp_lin2] * sum_MM;
            fit[temp_lin][temp_lin2] = fit_temp[temp_lin2];
            fit_BDR[temp_lin][temp_lin2] = fit_temp_MM[temp_lin2];
        }
    }

    for (temp_lin = 0; temp_lin < samplePoints; temp_lin++)
    {
        for (temp_lin2 = 0; temp_lin2 < samplePoints; temp_lin2++)
        {
            fit_trans[temp_lin][temp_lin2] = fit[temp_lin2][temp_lin];
        }
    }

    linspace(0, (samplePoints - 1) * z_min, samplePoints, temp_x);

    for (int temp_li = 0; temp_li < samplePoints; temp_li++)
    {
        int begin_inter, end_inter;
        double cos_temp;
        if (temp_li == 0)
        {
            begin_inter = 0;
            end_inter = c_photon_mask[temp_li];
        }
        else
        {
            begin_inter = c_photon_mask[temp_li - 1];
            end_inter = c_photon_mask[temp_li];
        }

        double S_temp[end_inter - begin_inter], sig_temp[end_inter - begin_inter];
        for (int temp_li2 = begin_inter; temp_li2 < end_inter; temp_li2++)
        {
            S_temp[temp_li2 - begin_inter] = DetS[temp_li2];
        }
        for (int temp_li3 = 0; temp_li3 < samplePoints; temp_li3++)
        {
            linear_interpolate(temp_x, fit_trans[temp_li3], samplePoints, S_temp, sig_temp, end_inter - begin_inter);
            for (int temp_li2 = begin_inter; temp_li2 < end_inter; temp_li2++)
            {
                if(!isnan(sqrt(sig_temp[temp_li2 - begin_inter] * DetW[temp_li2] * DetL[temp_li2])*cos_temp)){
                    cos_temp = cos(k_nolin[temp_li3] * S_temp[temp_li2 - begin_inter]);
                sig[temp_li3] += sqrt(sig_temp[temp_li2 - begin_inter] * DetW[temp_li2] * DetL[temp_li2]) * (cos_temp - I * sqrt(1.0 - cos_temp * cos_temp));
                }
            }
        }
    }
    printf("sig_name = %s\n", Miename);
    write_complex_array_to_binfile(Miename, sig, samplePoints);

    // /*single-scattered photons*/
    linspace(1, samplePoints, samplePoints, temp_x);

    for (temp_lin = 0; temp_lin < samplePoints; temp_lin++)
    {
        sum_fit = 0;
        sum_MM = 0;
        fit_polynomial(temp_x, MM_ss[temp_lin], samplePoints, coeffs);
        for (temp_lin2 = 0; temp_lin2 < samplePoints; temp_lin2++)
        {
            fit_temp[temp_lin2] = polyval(coeffs, ORDER, temp_x[temp_lin2]);
            if (fit_temp[temp_lin2] < 0)
                fit_temp[temp_lin2] = 0;
            sum_fit += fit_temp[temp_lin2];
            sum_MM += MM_ss[temp_lin][temp_lin2];
        }
        for (temp_lin2 = 0; temp_lin2 < samplePoints; temp_lin2++)
        {
            fit_temp[temp_lin2] = fit_temp[temp_lin2] / (sum_fit + ls);
            fit_temp_MM[temp_lin2] = fit_temp[temp_lin2] * sum_MM;
            fit[temp_lin][temp_lin2] = fit_temp[temp_lin2];
            fit_BDR[temp_lin][temp_lin2] = fit_temp_MM[temp_lin2];
        }
    }

    for (temp_lin = 0; temp_lin < samplePoints; temp_lin++)
    {
        for (temp_lin2 = 0; temp_lin2 < samplePoints; temp_lin2++)
        {
            fit_trans[temp_lin][temp_lin2] = fit[temp_lin2][temp_lin];
        }
    }

    linspace(0, (samplePoints - 1) * z_min, samplePoints, temp_x);

    for (int temp_li = 0; temp_li < samplePoints; temp_li++)
    {
        int begin_inter, end_inter;
        double cos_temp;
        if (temp_li == 0)
        {
            begin_inter = 0;
            end_inter = c_photon_mask_ss[temp_li];
        }
        else
        {
            begin_inter = c_photon_mask_ss[temp_li - 1];
            end_inter = c_photon_mask_ss[temp_li];
        }

        double S_temp[end_inter - begin_inter], sig_temp[end_inter - begin_inter];
        for (int temp_li2 = begin_inter; temp_li2 < end_inter; temp_li2++)
        {
            S_temp[temp_li2 - begin_inter] = DetS_ss[temp_li2];
        }
        for (int temp_li3 = 0; temp_li3 < samplePoints; temp_li3++)
        {
            linear_interpolate(temp_x, fit_trans[temp_li3], samplePoints, S_temp, sig_temp, end_inter - begin_inter);
            for (int temp_li2 = begin_inter; temp_li2 < end_inter; temp_li2++)
            {
                if(!isnan(sqrt(sig_temp[temp_li2 - begin_inter] * DetW_ss[temp_li2] * DetL_ss[temp_li2])*cos_temp)){
                    cos_temp = cos(k_nolin[temp_li3] * S_temp[temp_li2 - begin_inter]);
                sig_ss[temp_li3] += sqrt(sig_temp[temp_li2 - begin_inter] * DetW_ss[temp_li2] * DetL_ss[temp_li2]) * (cos_temp - I * sqrt(1.0 - cos_temp * cos_temp));
                }
            }
        }
    }
    printf("sig_ss_name = %s\n", Miename_ss);
    write_complex_array_to_binfile(Miename_ss, sig_ss, samplePoints);

    // /*multi-scattered photons*/
    linspace(1, samplePoints, samplePoints, temp_x);

    for (temp_lin = 0; temp_lin < samplePoints; temp_lin++)
    {
        sum_fit = 0;
        sum_MM = 0;
        fit_polynomial(temp_x, MM_ms[temp_lin], samplePoints, coeffs);
        for (temp_lin2 = 0; temp_lin2 < samplePoints; temp_lin2++)
        {
            fit_temp[temp_lin2] = polyval(coeffs, ORDER, temp_x[temp_lin2]);
            if (fit_temp[temp_lin2] < 0)
                fit_temp[temp_lin2] = 0;
            sum_fit += fit_temp[temp_lin2];
            sum_MM += MM_ms[temp_lin][temp_lin2];
        }
        for (temp_lin2 = 0; temp_lin2 < samplePoints; temp_lin2++)
        {
            fit_temp[temp_lin2] = fit_temp[temp_lin2] / (sum_fit + ls);
            fit_temp_MM[temp_lin2] = fit_temp[temp_lin2] * sum_MM;
            fit[temp_lin][temp_lin2] = fit_temp[temp_lin2];
            fit_BDR[temp_lin][temp_lin2] = fit_temp_MM[temp_lin2];
        }
    }

    for (temp_lin = 0; temp_lin < samplePoints; temp_lin++)
    {
        for (temp_lin2 = 0; temp_lin2 < samplePoints; temp_lin2++)
        {
            fit_trans[temp_lin][temp_lin2] = fit[temp_lin2][temp_lin];
        }
    }

    linspace(0, (samplePoints - 1) * z_min, samplePoints, temp_x);

    for (int temp_li = 0; temp_li < samplePoints; temp_li++)
    {
        int begin_inter, end_inter;
        double cos_temp;
        if (temp_li == 0)
        {
            begin_inter = 0;
            end_inter = c_photon_mask_ms[temp_li];
        }
        else
        {
            begin_inter = c_photon_mask_ms[temp_li - 1];
            end_inter = c_photon_mask_ms[temp_li];
        }

        double S_temp[end_inter - begin_inter], sig_temp[end_inter - begin_inter];
        for (int temp_li2 = begin_inter; temp_li2 < end_inter; temp_li2++)
        {
            S_temp[temp_li2 - begin_inter] = DetS_ms[temp_li2];
        }
        for (int temp_li3 = 0; temp_li3 < samplePoints; temp_li3++)
        {
            linear_interpolate(temp_x, fit_trans[temp_li3], samplePoints, S_temp, sig_temp, end_inter - begin_inter);
            for (int temp_li2 = begin_inter; temp_li2 < end_inter; temp_li2++)
            {
                if(!isnan(sqrt(sig_temp[temp_li2 - begin_inter] * DetW_ms[temp_li2] * DetL_ms[temp_li2])*cos_temp)){
                    cos_temp = cos(k_nolin[temp_li3] * S_temp[temp_li2 - begin_inter]);
                sig_ms[temp_li3] += sqrt(sig_temp[temp_li2 - begin_inter] * DetW_ms[temp_li2] * DetL_ms[temp_li2]) * (cos_temp - I * sqrt(1.0 - cos_temp * cos_temp));
                }
            }
        }
    }
    printf("sig_ms_name = %s\n", Miename_ms);
    write_complex_array_to_binfile(Miename_ms, sig_ms, samplePoints);
    

    printf("Postprocessing complete.\n");
    now = time(NULL);
    printf("%s\n", ctime(&now));
    free(DetS);
    free(DetL);
    free(DetW);
    free(DetS_ss);
    free(DetL_ss);
    free(DetW_ss);
    free(DetS_ms);
    free(DetL_ms);
    free(DetW_ms);
    for (int i = 0; i < particle_center_num; i++) {
        free(particles_center[i]);
    }
    free(particles_center);
    free(particle_center_index);
}
    for (int i = 0; i < particle_num; i++) {
        free(particles[i]);
    }
    free(particles);
    free(particle_index);
    return 0;
} /* end of main */

/* SUBROUTINES */

/**************************************************************************
 *	RandomGen
 *      A random number generator that generates uniformly
 *      distributed random numbers between 0 and 1 inclusive.
 *      The algorithm is based on:
 *      W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P.
 *      Flannery, "Numerical Recipes in C," Cambridge University
 *      Press, 2nd edition, (1992).
 *      and
 *      D.E. Knuth, "Seminumerical Algorithms," 2nd edition, vol. 2
 *      of "The Art of Computer Programming", Addison-Wesley, (1981).
 *
 *      When Type is 0, sets Seed as the seed. Make sure 0<Seed<32000.
 *      When Type is 1, returns a random number.
 *      When Type is 2, gets the status of the generator.
 *      When Type is 3, restores the status of the generator.
 *
 *      The status of the generator is represented by Status[0..56].
 *
 *      Make sure you initialize the seed before you get random
 *      numbers.
 ****/
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC 1.0E-9

double RandomGen(char Type, long Seed, long *Status)
{
    static long i1, i2, ma[56]; /* ma[0] is not used. */
    long mj, mk;
    short i, ii;

    if (Type == 0)
    { /* set seed. */
        mj = MSEED - (Seed < 0 ? -Seed : Seed);
        mj %= MBIG;
        ma[55] = mj;
        mk = 1;
        for (i = 1; i <= 54; i++)
        {
            ii = (21 * i) % 55;
            ma[ii] = mk;
            mk = mj - mk;
            if (mk < MZ)
                mk += MBIG;
            mj = ma[ii];
        }
        for (ii = 1; ii <= 4; ii++)
            for (i = 1; i <= 55; i++)
            {
                ma[i] -= ma[1 + (i + 30) % 55];
                if (ma[i] < MZ)
                    ma[i] += MBIG;
            }
        i1 = 0;
        i2 = 31;
    }
    else if (Type == 1)
    { /* get a number. */
        if (++i1 == 56)
            i1 = 1;
        if (++i2 == 56)
            i2 = 1;
        mj = ma[i1] - ma[i2];
        if (mj < MZ)
            mj += MBIG;
        ma[i1] = mj;
        return (mj * FAC);
    }
    else if (Type == 2)
    { /* get status. */
        for (i = 0; i < 55; i++)
            Status[i] = ma[i + 1];
        Status[55] = i1;
        Status[56] = i2;
    }
    else if (Type == 3)
    { /* restore status. */
        for (i = 0; i < 55; i++)
            ma[i + 1] = Status[i];
        i1 = Status[55];
        i2 = Status[56];
    }
    else
        puts("Wrong parameter to RandomGen().");
    return (0);
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

/***********************************************************
 *  Determine if the two position are located in the same voxel
 *	Returns 1 if same voxel, 0 if not same voxel.
 ****/
Boolean SameVoxel(double x1, double y1, double z1, double x2, double y2, double z2,
                  double dx, double dy, double dz)
{
    double xmin = min2((floor)(x1 / dx), (floor)(x2 / dx)) * dx;
    double ymin = min2((floor)(y1 / dy), (floor)(y2 / dy)) * dy;
    double zmin = min2((floor)(z1 / dz), (floor)(z2 / dz)) * dz;
    double xmax = xmin + dx;
    double ymax = ymin + dy;
    double zmax = zmin + dz;
    Boolean sv = 0;

    sv = (x1 <= xmax && x2 <= xmax && y1 <= ymax && y2 <= ymax && z1 < zmax && z2 <= zmax);
    return (sv);
}

/***********************************************************
 * max2
 ****/
double max2(double a, double b)
{
    double m;
    if (a > b)
        m = a;
    else
        m = b;
    return m;
}

/***********************************************************
 * min2
 ****/
double min2(double a, double b)
{
    double m;
    if (a >= b)
        m = b;
    else
        m = a;
    return m;
}
/***********************************************************
 * min3
 ****/
double min3(double a, double b, double c)
{
    double m;
    if (a <= min2(b, c))
        m = a;
    else if (b <= min2(a, c))
        m = b;
    else
        m = c;
    return m;
}

/* How much step size will the photon take to get the first voxel crossing in one single
    long step? */
//
double FindVoxelFace2(int *fdir, double x1, double y1, double z1, int *det_num, int Pick_det, double detx, double det_radius, double det_z, double cos_accept, int Ndetectors, double dx, double dy, double dz, double ux, double uy, double uz)
{

    int ix1 = floor(x1 / dx);
    int iy1 = floor(y1 / dy);
    int iz1 = floor(z1 / dz);
    int izd = floor(det_z / dz);

    // ix2, iy2, iz2: indices of the voxel faces lying ahead of the photon's propagation path
    int ix2, iy2, iz2;
    // equal to equation 4.12 in Zhao's thesis
    if (ux >= 0)
        ix2 = ix1 + 1;
    else
        ix2 = ix1;
    // equal to equation 4.13 in Zhao's thesis
    if (uy >= 0)
        iy2 = iy1 + 1;
    else
        iy2 = iy1;
    // equal to equation 4.14 in Zhao's thesis
    if (uz >= 0)
        iz2 = iz1 + 1;
    else
        iz2 = iz1;
    // xs, ys, zs: distances from these voxel faces to the current position of the photon utilizing its propagation directions
    double xs = fabs((ix2 * dx - x1) / ux);
    double ys = fabs((iy2 * dy - y1) / uy);
    double zs = fabs((iz2 * dz - z1) / uz);
    // s: desired distance of the photon to its closest voxel face
    double s = min3(xs, ys, zs);
    *fdir = s == zs? 3: s == ys? 2: 1;
    // check detection
    if (-uz >= cos_accept && izd == iz1 && s == zs && fabs(y1 + s * uy) <= det_radius && fabs(-uz - cos_accept)<=2e-4)
    {
        if (fabs(x1 + s * ux - detx) <= det_radius)
            *det_num = Pick_det;
    }
    return (s);
}

double TT_final1(double gf, double gb, double alpha_f, double alpha_b, double C, double ang)
{
    double kf = alpha_f * gf / PI * pow(1-gf,2*alpha_f) * pow(1+gf,2*alpha_f) / (pow(1+gf,2*alpha_f) - pow(1-gf,2*alpha_f));
    double kb = alpha_b * gb / PI * pow(1-gb,2*alpha_b) * pow(1+gb,2*alpha_b) / (pow(1+gb,2*alpha_b) - pow(1-gb,2*alpha_b));
    double core_f = kf*pow((1 + gf * gf -2 * gf * ang), -alpha_f-1);
    double core_b = kb*pow((1 + gb * gb -2 * gb * (-ang)), -alpha_b-1);
    double TT_final = C * core_f + (1-C) * core_b;

    return TT_final;
}

double TT_is_final1(double gf, double gb, double alpha_f, double alpha_b, double C, double ang)
{
    double kf = alpha_f * gf / PI * pow(1-gf,2*alpha_f) * pow(1+gf*gf,alpha_f) / (pow(1+gf*gf,alpha_f) - pow(1-gf,2*alpha_f));
    double kb = alpha_b * gb / PI * pow(1-gb,2*alpha_b) * pow(1+gb*gb,alpha_b) / (pow(1+gb*gb,alpha_b) - pow(1-gb,2*alpha_b));
    double core_f = kf*pow((1 + gf * gf -2 * gf * ang), -alpha_f-1);
    double core_b = kb*pow((1 + gb * gb -2 * gb * (1-ang)), -alpha_b-1);
    double TT_final = C * core_f + (1-C) * core_b;

    return TT_final;
}

double polyval(double *coeffs, int order, double x)
{
    double result = 0.0;
    int i;
    for (i = order; i >= 0; i--)
    {
        result = result * x + coeffs[i];
    }
    return result;
}

void fit_polynomial(double *x, double *y, int num_points, double *coeffs)
{
    double A[ORDER + 1][ORDER + 1] = {0};
    double b[ORDER + 1] = {0};
    int i, j, k;

    // construct the A matrix and b vector
    for (i = 0; i < num_points; i++)
    {
        for (j = 0; j <= ORDER; j++)
        {
            for (k = 0; k <= ORDER; k++)
            {
                A[j][k] += pow(x[i], j + k);
            }
            b[j] += pow(x[i], j) * y[i];
        }
    }

    // solve for the coefficients using Gaussian elimination
    for (k = 0; k < ORDER; k++)
    {
        for (i = k + 1; i <= ORDER; i++)
        {
            double factor = A[i][k] / (A[k][k]+ls);
            for (j = k; j <= ORDER; j++)
            {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }
    }
    for (i = ORDER; i >= 0; i--)
    {
        coeffs[i] = b[i];
        for (j = i + 1; j <= ORDER; j++)
        {
            coeffs[i] -= A[i][j] * coeffs[j];
        }
        coeffs[i] /= A[i][i];
    }
}

void linspace(double start, double stop, int num, double *arr)
{
    double step = (stop - start) / (double)(num - 1);
    for (int i = 0; i < num; i++)
    {
        arr[i] = start + i * step;
    }
}

void linear_interpolate(double *x_array, double *y_array, int size_x, double *x_interp, double *y_interp, int size_interp)
{
    int i, j;
    double x1, y1, x2, y2, slope, y_intercept;

    // Iterate over the interpolation points
    for (i = 0; i < size_interp; i++)
    {
        // Find the index of the left and right points for the interpolation point
        for (j = 0; j < size_x - 1; j++)
        {
            if (x_interp[i] < x_array[j + 1])
            {
                break;
            }
        }
        // Use linear interpolation to find the y-value at the interpolation point
        x1 = x_array[j];
        y1 = y_array[j];
        x2 = x_array[j + 1];
        y2 = y_array[j + 1];
        slope = (y2 - y1) / (x2 - x1+ls);
        y_intercept = y1 - slope * x1;
        y_interp[i] = slope * x_interp[i] + y_intercept;
    }
}

void write_complex_array_to_binfile(char *filename, double complex *data, int size)
{
    FILE *f = fopen(filename, "wb");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    fwrite(data, sizeof(double complex), size, f);

    fclose(f);
}

double distanceToLineSegment(Point3D point, Point3D start, Point3D end) {
    double dirX = end.x - start.x;
    double dirY = end.y - start.y;
    double dirZ = end.z - start.z;

    double vecX = point.x - start.x;
    double vecY = point.y - start.y;
    double vecZ = point.z - start.z;

    double vecEndX = point.x - end.x;
    double vecEndY = point.y - end.y;
    double vecEndZ = point.z - end.z;

    double distanceStart = sqrt(vecX * vecX + vecY * vecY + vecZ * vecZ);
    double distanceEnd = sqrt(vecEndX * vecEndX + vecEndY * vecEndY + vecEndZ * vecEndZ);

    double projection = (vecX * dirX + vecY * dirY + vecZ * dirZ) / (dirX * dirX + dirY * dirY + dirZ * dirZ);

    if (projection < 0.0) {
        return distanceStart;
    }
    else if (projection > 1.0) {
        return distanceEnd;
    }
    else {
        double projX = start.x + projection * dirX;
        double projY = start.y + projection * dirY;
        double projZ = start.z + projection * dirZ;

        double distanceProjection = sqrt((point.x - projX) * (point.x - projX) +
                                         (point.y - projY) * (point.y - projY) +
                                         (point.z - projZ) * (point.z - projZ));

        return distanceProjection;
    }
}

void reshapeArray(const double* input, int rows, int cols, double*** output) {
    *output = (double**)malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        (*output)[i] = (double*)malloc(cols * sizeof(double));
        for (int j = 0; j < cols; j++) {
            (*output)[i][j] = input[i * cols + j];
        }
    }
}

void findIndices(const long* array, int size, long target, int* lowerIndex, int* upperIndex) {
    int left = 0;
    int right = size - 1;

    while (left <= right) {
        int mid = left + (right - left) / 2;

        if (array[mid] <= target) {
            *lowerIndex = mid;
            left = mid + 1;
        } else {
            *upperIndex = mid;
            right = mid - 1;
        }
    }
}

double RFresnel(double n1,		/* incident refractive index.*/
                double n2,		/* transmit refractive index.*/
                double ca1,		/* cosine of the incident */
                /* angle a1, 0<a1<90 degrees. */
                double *ca2_Ptr) 	/* pointer to the cosine */
/* of the transmission */
/* angle a2, a2>0. */
{
    double r;
    
    if(n1==n2) { /** matched boundary. **/
        *ca2_Ptr = ca1;
        r = 0.0;
	}
    else if(ca1>(1.0 - 1.0e-12)) { /** normal incidence. **/
        *ca2_Ptr = ca1;
        r = (n2-n1)/(n2+n1);
        r *= r;
	}
    else if(ca1< 1.0e-6)  {	/** very slanted. **/
        *ca2_Ptr = 0.0;
        r = 1.0;
	}
    else  {			  		/** general. **/
        double sa1, sa2; /* sine of incident and transmission angles. */
        double ca2;      /* cosine of transmission angle. */
        sa1 = sqrt(1-ca1*ca1);
        sa2 = n1*sa1/n2;
        if(sa2>=1.0) {	
            /* double check for total internal reflection. */
            *ca2_Ptr = 0.0;
            r = 1.0;
		}
        else {
            double cap, cam;	/* cosines of sum ap or diff am of the two */
            /* angles: ap = a1 + a2, am = a1 - a2. */
            double sap, sam;	/* sines. */
            *ca2_Ptr = ca2 = sqrt(1-sa2*sa2);
            cap = ca1*ca2 - sa1*sa2; /* c+ = cc - ss. */
            cam = ca1*ca2 + sa1*sa2; /* c- = cc + ss. */
            sap = sa1*ca2 + ca1*sa2; /* s+ = sc + cs. */
            sam = sa1*ca2 - ca1*sa2; /* s- = sc - cs. */
            r = 0.5*sam*sam*(cam*cam+cap*cap)/(sap*sap*cam*cam); 
            /* rearranged for speed. */
		}
	}
    return(r);
} /******** END SUBROUTINE **********/


void CrossUpOrNot(double ni, double nt){
                            
    double cos_ctheta = ni>nt ? sqrt(1.0 - nt*nt/(ni*ni)) : 0.0;
    double r = 0.0;
    double uz1;
    if(-uz <= cos_ctheta)
        r = 1.0;
    else
        r = RFresnel(ni, nt, -uz, &uz1);
    /*do transmission*/
    rnd = RandomNum;

    
	if(rnd > r){
        if(iz==0)  {
            uz = -uz1;
            photon_status = DEAD;
        }
        else {
            ux *= ni/nt;
            uy *= ni/nt;
            uz = -uz1;
        }
    }
    
    else{
        iz += 1;
        // uz = -uz;
        sleft += 2*ls;
        calculateReflection(&ux, &uy, &uz, 
                         ux, uy, uz, 
                         x, y, z, 
                         0, 0,0.08); 
        haven_scatter=1;
        zfocus_used = 2*z-zfocus_used;
    }
        
}

void CrossDnOrNot(double ni, double nt){
                            
    double cos_ctheta = ni>nt ? sqrt(1.0 - nt*nt/(ni*ni)) : 0.0;
    double r = 0.0;
    double uz1;
    if(uz <= cos_ctheta)
        r = 1.0;
    else
        r = RFresnel(ni, nt, uz, &uz1);
    /*do transmission*/
    // while ((rnd = RandomNum) <= 0.0);   /* yields 0 < rnd <= 1 */
    rnd = RandomNum;
	if(rnd > r){
        if(iz == Nz - 1)  {
            uz = uz1;
            photon_status = DEAD;
        }
        else {
            ux *= ni/nt;
            uy *= ni/nt;
            uz = uz1;
        }
    }
    
    else{
        iz -= 1;
        // uz = -uz;
        sleft += 2*ls;
        zfocus_used = 2*z-zfocus_used;
        calculateReflection(&ux, &uy, &uz, 
                         ux, uy, uz, 
                         x, y, z, 
                         0, 0,0.08); 
        haven_scatter=1;
    }
        
}

void CrossOutOrNot(double ni, double nt){
                            
    double cos_ctheta = ni>nt ? sqrt(1.0 - nt*nt/(ni*ni)) : 0.0;
    double r = 0.0;
    double uy1;
    if(-uy <= cos_ctheta)
        r = 1.0;
    else
        r = RFresnel(ni, nt, -uy, &uy1);
    /*do transmission*/
    rnd = RandomNum;
    
	if(rnd > r){
        if(iy==0)  {
            uy = -uy1;
            photon_status = DEAD;
        }
        else {
            ux *= ni/nt;
            uz *= ni/nt;
            uy = -uy1;
        }
    }
    
    else{
        iy += 1;
        // uy = -uy;
        calculateReflection(&ux, &uy, &uz, 
                         ux, uy, uz, 
                         x, y, z, 
                         0, 0,0.08); 
        sleft += 2*ls;
        haven_scatter=1;
    }
        
}

void CrossInOrNot(double ni, double nt){
                            
    double cos_ctheta = ni>nt ? sqrt(1.0 - nt*nt/(ni*ni)) : 0.0;
    double r = 0.0;
    double uy1;
    if(uy <= cos_ctheta)
        r = 1.0;
    else
        r = RFresnel(ni, nt, uy, &uy1);
    /*do transmission*/
    rnd = RandomNum;
    
	if(rnd > r){
        if(iy == Ny - 1)  {
            uy = uy1;
            photon_status = DEAD;
        }
        else {
            ux *= ni/nt;
            uz *= ni/nt;
            uy = uy1;
        }
    }
    
    else{
        iy -= 1;
        // uy = -uy;
        calculateReflection(&ux, &uy, &uz, 
                         ux, uy, uz, 
                         x, y, z, 
                         0, 0,0.08); 
        sleft += 2*ls;
        haven_scatter=1;
    }
        
}

void CrossLeftOrNot(double ni, double nt){
                            
    double cos_ctheta = ni>nt ? sqrt(1.0 - nt*nt/(ni*ni)) : 0.0;
    double r = 0.0;
    double ux1;
    if(-ux <= cos_ctheta)
        r = 1.0;
    else
        r = RFresnel(ni, nt, -ux, &ux1);
    /*do transmission*/
    rnd = RandomNum;
    
	if(rnd > r){
        if(ix==0)  {
            ux = -ux1;
            photon_status = DEAD;
        }
        else {
            uz *= ni/nt;
            uy *= ni/nt;
            ux = -ux1;
        }
    }
    
    else{
        ix += 1;
        // ux = -ux;
        calculateReflection(&ux, &uy, &uz, 
                         ux, uy, uz, 
                         x, y, z, 
                         0, 0,0.08); 
        sleft += 2*ls;
        haven_scatter=1;
    }
        
}

void CrossRightOrNot(double ni, double nt){
                            
    double cos_ctheta = ni>nt ? sqrt(1.0 - nt*nt/(ni*ni)) : 0.0;
    double r = 0.0;
    double ux1;
    if(ux <= cos_ctheta)
        r = 1.0;
    else
        r = RFresnel(ni, nt, ux, &ux1);
    /*do transmission*/
    rnd = RandomNum;
    
	if(rnd > r){
        if(ix == Nx - 1)  {
            ux = ux1;
            photon_status = DEAD;
        }
        else {
            uz *= ni/nt;
            uy *= ni/nt;
            ux = ux1;
        }
    }
    
    else{
        ix -= 1;
        // ux = -ux;
        calculateReflection(&ux, &uy, &uz, 
                         ux, uy, uz, 
                         x, y, z, 
                         0, 0,0.08); 
        sleft += 2*ls;
        haven_scatter=1;
    }
        
}
void generateRandomGaussian(double *x, double *y, double r1, double r2) {
    double u1, u2;
    double r, theta;

    u1 = r1;
    u2 = r2;

    r = sqrt(-1.0 * log(u1));
    theta = 2.0 * PI * u2;

    *x = r * cos(theta);
    *y = r * sin(theta);
}

void generateGaussianDirection(double x, double y, double z, double *ux, double *uy, double *uz, double zr, double zfocus, double detx) {
    double R_z;
    double delta_z = z - zfocus;
    x = x-detx;
    R_z = -(delta_z)*sqrt(1+(zr/(delta_z))*(zr/(delta_z)));
    if(R_z==0){
        *ux = 0;
        *uy = 0;
        *uz = 1.0;
    }
    else{
        double d_amp = 1.0 / sqrt(1 + (x * x + y * y) / (R_z * R_z));
        *ux = -x / R_z * d_amp;
        *uy = -y / R_z * d_amp;
        *uz = 1.0 * d_amp;
    }
    
}

void calculateReflection(double* reflectedX, double* reflectedY, double* reflectedZ, 
                         double incidentX, double incidentY, double incidentZ, 
                         double pointAX, double pointAY, double pointAZ, 
                         double pointBX, double pointBY, double pointBZ) {
    double normalX = pointBX - pointAX;
    double normalY = pointBY - pointAY;
    double normalZ = pointBZ - pointAZ;

    double length = sqrt(normalX * normalX + normalY * normalY + normalZ * normalZ);
    normalX /= length;
    normalY /= length;
    normalZ /= length;

    length = sqrt(incidentX * incidentX + incidentY * incidentY + incidentZ * incidentZ);
    incidentX /= length;
    incidentY /= length;
    incidentZ /= length;

    double dotProduct = incidentX * normalX + incidentY * normalY + incidentZ * normalZ;
    *reflectedX = incidentX - 2 * dotProduct * normalX;
    *reflectedY = incidentY - 2 * dotProduct * normalY;
    *reflectedZ = incidentZ - 2 * dotProduct * normalZ;

    length = sqrt(*reflectedX * *reflectedX + *reflectedY * *reflectedY + *reflectedZ * *reflectedZ);
    *reflectedX /= length;
    *reflectedY /= length;
    *reflectedZ /= length;
}

double cos_final1(double gf, double gb, double alpha_f, double alpha_b, double C)
{
    double q1f = (1+gf*gf);
    double q2f = powf(1-gf,2*alpha_f);
    double q3f = powf(1+gf,2*alpha_f);
    double q1b = (1+gb*gb);
    double q2b = powf(1-gb,2*alpha_b);
    double q3b = powf(1+gb,2*alpha_b);

    rnd = RandomNum;
    if (rnd<=C){ /*more common event: forward scatter (C close to 0.9) */
    rnd = RandomNum;
    costheta = (q1f/(2*gf))-powf(rnd/q2f + (1-rnd)/q3f,-1/alpha_f)/(2*gf);
    }
    else{ /* more rare event: backward scatter */
    rnd = RandomNum;
    costheta = (q1b/(2*gb))-powf(rnd/q2b + (1-rnd)/q3b,-1/alpha_b)/(2*gb);
    costheta = -costheta;
    }
    if (costheta>=1.0) costheta = 1-1e-9;
    if (costheta<=-1.0) costheta = -(1-1e-9);
    

    return costheta;
}

double cos_is_final1(double gf, double gb, double alpha_f, double alpha_b, double C)
{
    double q1f = (1+gf*gf);
    double q2f = powf(1-gf,2*alpha_f);
    double q3f = powf(1+gf*gf,alpha_f);
    double q1b = (1+gb*gb);
    double q2b = powf(1-gb,2*alpha_b);
    double q3b = powf(1+gb*gb,alpha_b);

    rnd = RandomNum;
    if (rnd<=C){ /*more common event: forward scatter (C close to 0.9) */
    rnd = RandomNum;
    costheta = (q1f/(2*gf))-powf(rnd/q2f + (1-rnd)/q3f,-1/alpha_f)/(2*gf);
    }
    else{ /* more rare event: backward scatter */
    rnd = RandomNum;
    costheta = (q1b/(2*gb))-powf(rnd/q2b + (1-rnd)/q3b,-1/alpha_b)/(2*gb);
    costheta = 1-costheta;
    }
    if (costheta>=1.0) costheta = 1-1e-9;
    if (costheta<=-1.0) costheta = -(1-1e-9);
    

    return costheta;
}