### Official implementation of "Monte Carlo-based realistic simulation of optical coherence tomography angiography"

# Example
If you want to obtain similar results like those in our paper, please download the examples [here](https://drive.google.com/file/d/17Pwpd4Wvzu0sZblR-NfzYP1sEmP-644f/view?usp=drive_link), and run the lookImage.m in each dictionary.

<img src="https://github.com/Jianing-Mao/OCTA_MC/blob/master/example/Bscan.png" width="500px">
<img src="https://github.com/Jianing-Mao/OCTA_MC/blob/master/example/Bscan2.png" width="500px">
<img src="https://github.com/Jianing-Mao/OCTA_MC/blob/master/example/qocta.png" width="500px">

## Get Started
The code is builded in MATLAB R2020a and ANSI C. Please install the [MatScat](https://ww2.mathworks.cn/matlabcentral/fileexchange/36831-matscat) in your MATLAB.

## Begin Monte Carlo simulation
### Step I: prepare the parameters

Run [TT_para_gen.m](https://github.com/Jianing-Mao/OCTA_MC/blob/master/TT_para_gen.m) to generate the wavelength-dependent parameters for TT scattering phase functions.

After this, you can get text files: _parameter_#.txt_ in dictionary [parameters](https://github.com/Jianing-Mao/OCTA_MC/tree/master/parameters), # is the diameter of the particle.

Run [mus_gen.m](https://github.com/Jianing-Mao/OCTA_MC/blob/master/mus_gen.m) to generate the wavelength-dependent scattering coefficients and efficient particle radii for static and dynamic media.

After this, you can get .mat files: _us\_#\_1300.mat_ and _r\_#.mat_ in dictionary [parameters](https://github.com/Jianing-Mao/OCTA_MC/tree/master/parameters), # is the diameter of the particle.

Run [GenTissue.m](https://github.com/Jianing-Mao/OCTA_MC/blob/master/GenTissue.m) to generate the simulation setting at all wavelengths.

After this, you can get _M_ _H.mci files and one _T.bin file in dictionary [settings](https://github.com/Jianing-Mao/OCTA_MC/tree/master/settings).

Run [Gen_particles.m](https://github.com/Jianing-Mao/OCTA_MC/blob/master/Gen_particles.m) to generate the random positions of the particles and the flows.

After this, you can get N+1 .bin file: _data_static_particles.bin_, _data_dynamic_particles_init.bin_, _data_dynamic_particles\_*.bin_ in dictionary [particles](https://github.com/Jianing-Mao/OCTA_MC/tree/master/particles), * in [2, N], N is the number of bscans for OCTA;

### Step II: run the Monte Carlo

Compile the [fullwaveMC_Bscan_OCTA.c](https://github.com/Jianing-Mao/OCTA_MC/blob/master/fullwaveMC_Bscan_OCTA.c)
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

Run [lookImage.m](https://github.com/Jianing-Mao/OCTA_MC/blob/master/lookImage.m) to show the OCT, OCTA and quantitative OCTA results.

**Make sure the parameters (wavelength, diameter, RI, number density, Bscan_number...) between the "Params settings START" and "Params settings DONE" are the same in every .m files.**

# To be implemented
* Integrate all codes into a calling function (now we need to carefully check the parameters of the generation, running, and viewing codes to make sure that all the parameters are the same.)
* Heterogeneous medium
* Mesh-based phantom model
* ...

# Acknowledgement
The codes are built on [MatScat](https://ww2.mathworks.cn/matlabcentral/fileexchange/36831-matscat), open-source codes in [OMLC](https://omlc.org/software/mc/) and [OCT_MC](https://github.com/RMTariant/OCT_MC). We sincerely appreciate the authors for sharing their codes. If you also want to try a more classical (such as using pencil beams and random sampling step sizes) but still a full-spectrum simulation, please refer to [fullwaveOCT code](https://github.com/Jianing-Mao/fullwaveOCT) and [paper](https://opg.optica.org/boe/fulltext.cfm?uri=boe-14-9-4644&id=536404)

# Contact
If you have any questions or any comments about the simulation, please do not hesitate to contact [jianing.mao@sjtu.edu.cn](jianing.mao@sjtu.edu.cn)
