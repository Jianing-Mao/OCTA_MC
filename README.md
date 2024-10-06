### Official implementation of "Monte Carlo-based realistic simulation of optical coherence tomography angiography"

## Get Started
The code is builded in MATLAB R2020a and ANSI C. Please install the [MatScat](https://ww2.mathworks.cn/matlabcentral/fileexchange/36831-matscat) in your MATLAB.

## Prepare scattering parameters
Run [TT_para_gen.m](https://github.com/Jianing-Mao/OCTA_MC/blob/master/Para_gen/TT_para_gen.m) to generate the wavelength-dependent parameters for TT scattering phase functions.

After this, you can get a text file: _parameter_#.txt_, # is the diameter of the particle; you should run multiple times to generate parameters for different types of media.

Run [mus_gen.m](https://github.com/Jianing-Mao/OCTA_MC/blob/master/Para_gen/mus_gen.m) to generate the wavelength-dependent scattering coefficients and efficient particle radii for static and dynamic media.

After this, you can get two .mat file: _us\_#\_1300.mat_ and _r\_#.mat_, # is the diameter of the particle; you should run multiple times to generate parameters for different types of media.

Run [Gen_particles.m](https://github.com/Jianing-Mao/OCTA_MC/blob/master/Para_gen/Gen_particles.m) to generate the random positions of the particles and the flows.

After this, you can get N+1 .bin file: _data_static_particles.bin_, _data_dynamic_particles_init.bin_, _data_dynamic_particles\_*.bin_, * in [2, N], N is the number of bscans for OCTA;

Make sure the parameters (wavelength, angles, diameter, RI, number density, Bscan_number...) are the same in every files. Copy all the outputs to the [code](https://github.com/Jianing-Mao/OCTA_MC/tree/master/Code) directory.
## Begin Monte Carlo simulation
### Step I: prepare the tissue

Run [GenTissue.m](https://github.com/Jianing-Mao/OCTA_MC/blob/master/Code/GenTissue.m) to generate the simulation setting at all wavelengths.

### Step II: run the Monte Carlo

Compile the [fullwaveMC_Bscan_OCTA.c](https://github.com/Jianing-Mao/OCTA_MC/blob/master/Code/fullwaveMC_Bscan_OCTA.c)
```sh
gcc fullwaveMC_Bscan_OCTA.c -lm -o test
```

Run the generated file with two arguments: (1) arg1: the index of the detectors; (2) arg2: the number of Bscans:
```sh
./test 256 4
```
You can generate a Bscan by running the simulation of all of the detector in parallel:
```sh
./test 1 4
./test 2 4
...
```
### Step III: show the results

Run [lookImage.m](https://github.com/Jianing-Mao/OCTA_MC/blob/master/Code/lookImage.m) to show the OCT, OCTA and quantitative OCTA results.

# To be implemented
* Integrate all parameters and procedures into a calling function (now we need to carefully check the parameters of the generation, running, and viewing codes to make sure that all the parameters are the same.)
* Heterogeneous medium
* Mesh-based phantom model
* ...

# Example
<img src="https://github.com/Jianing-Mao/OCTA_MC/blob/master/example/Bscan.png" width="500px">
<img src="https://github.com/Jianing-Mao/OCTA_MC/blob/master/example/qocta.png" width="500px">

# Acknowledgement
The codes are built on [MatScat](https://ww2.mathworks.cn/matlabcentral/fileexchange/36831-matscat), open-source codes in [OMLC](https://omlc.org/software/mc/) and [OCT_MC](https://github.com/RMTariant/OCT_MC). We sincerely appreciate the authors for sharing their codes. If you also want to try a more classical (such as using pencil beams and random sampling step sizes) but still a full-spectrum simulation, please refer to [fullwaveOCT code](https://github.com/Jianing-Mao/fullwaveOCT) and [paper](https://opg.optica.org/boe/fulltext.cfm?uri=boe-14-9-4644&id=536404)
