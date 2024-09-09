# TSE

## Triple stellar evolution code

Single stellar evolution and binary interactions are modelled using the binary population synthesis code MOBSE (https://mobse-webpage.netlify.app) which builds upon the code BSE (https://www.ascl.net/1303.014). We have slighlty modified the files `mobse/src/evolve.f`, `mobse/src/mobse.f`, and `mobse/input/const_mobse.h` in MOBSE. These changes are provided as a patch (see below) and were only made to ensure compatibility between TSE and MOBSE, e.g., change of array lengths and MOBSE input/output, but the physical prescription of MOBSE remains unaltered.  Using other versions of MOBSE or BSE-type codes is in general possible, as long as both modified files are maintained.

TSE evolves the secular dynamics of hierarchical triples and combines it with the stellar evolution of the inner binary and the tertiary companion.

### Installation

In order to install the code follow these steps to install

1. Add MOBSE, e.g., by cloning `git clone git@gitlab.com:mobse/source-code.git mobse` or by going to https://mobse-webpage.netlify.app and downloading the MOBSE directory. In any case, MOBSE should be added to the root directory be named `mobse`.
2. Navigate to the `mobse` directory: `cd mobse`.
3. Implement our changes to `mobse` by `patch -p1 < ../mobse.patch`. This should notify you about patching the files `mobse/src/evolve.f`, `mobse/src/mobse.f`, and `mobse/input/const_mobse.h`.
4. Navigate to the `mobse/src` directory: `cd src`.
5. Execute the `Makefile` by running `make mobse`. You may need to change `mobse/src/Makefile` according to the specifications of your fortran compiler.
6. Navigate back to the root directory: `cd ./../..`.
7. Create two directories where the outcome of the simulation is stored, `mkdir plots` and `mkdir data`.
8. Install all required packages by executing `pip install -r requirements.txt`.

Alternatively, simply execute `compile.sh` in the root directory. You may want to change the `mobse/src/Makefile` according to the specifications of your fortran compiler.

### Running

To run a system activate the environment (see above) and go to the `bin` directory. Then execute

`python tse.py`.

You may want to check out the number of optional arguments that are available by running

`python tse.py -h`.

E.g., do 

`python tse.py --random TRUE --seed 1 --Z 0.02`

to evolve a triple from a random distribution at solar metallicity.

If you are looking whether a triple leads to a particular configuration (e.g., the formation of BH triples) check out the function `CustomEvent` in `src/main.py` (including a working example).

TSE runs stably with the python version and packages listed in `requirements.txt`.

### Testing and output

Let's test the code and inspect the output by running, e.g., `python tse.py --random TRUE --seed 139 --Z 0.0001 --stellar_tides TRUE --max_time 15000.0 --method DOP853 --bhflag 1 --nsflag 3 --lamb 0.1 --alpha1 0.2`. This will lead to a terminal output similar as following.

First, the sampled initial conditions are printed:

```
Initial conditions:
 Seed=139,
 Z=0.000100,
 m1=34.868610,
 m2=31.640159,
 m3=24.505861,
 a_in=4.320278e+04,
 e_in=0.041194,
 a_out=1.726712e+06,
 e_out=0.759993,
 cos_i_in=-0.668036,
 omega_in=1.832684,
 Omega_in=3.141593,
 cos_i_out=-0.164133,
 omega_out=6.104407,
 Omega_out=0.000000
```

Then, we get a standard MOBSE/BSE output tabulating how the inner binary would have evolved if there was no tertiary companion:

```
Evolve inner binary as if it was isolated


      TIME      M1       M2   K1 K2        SEP    ECC  R1/ROL1 R2/ROL2  TYPE
     0.0000   34.86861   31.64016  1  1    43187.60156  0.04   0.000   0.000  INITIAL 
     5.8129   34.83288   31.61894  2  1    43224.61719  0.04   0.001   0.001  KW_CHNGE
     5.8199   34.83244   31.61883  4  1    43224.97266  0.04   0.003   0.001  KW_CHNGE
     6.2848   34.65763   31.61121  4  2    43343.74609  0.04   0.040   0.001  KW_CHNGE
     6.2930   34.63868   31.61083  4  4    43356.28125  0.04   0.046   0.002  KW_CHNGE
     6.3532   34.40531   31.60643  5  4    43508.64453  0.04   0.082   0.002  KW_CHNGE
     6.3644   33.87917   31.60555 14  4    43841.69531  0.04   0.000   0.002  KW_CHNGE
     6.7871   33.87919   31.50193 14  4    43886.07031  0.04   0.000   0.021  BEG_SYMB
     6.8721   33.87977   31.33843 14  5    43994.22266  0.04   0.000   0.076  KW_CHNGE
     6.8841   33.88026   29.58524 14 14    32407.31641  0.77   0.000   0.000  KW_CHNGE
 15000.0000   33.88026   29.58524 14 14    32388.83203  0.77   0.000   0.000  MAX_TIME
 Fallback:  0.95998213867552118
```

At next, the triple is evolved. For this purpose, the stellar evolution of each star, i.e., the evolution of its mass, radius, etc., is tabulated (first the primary, then the secondary, and finally the tertiary). This is achieved by evolving each of the stars with MOBSE with some artificial low-mass black hole companion at unphysical large distance, i.e., such that the stars are effectively single. Therefore, ignore the entries for M2, K2, SEP, ECC, R1/ROL1, R2/ROL2. All we care about are M1 (and other parameters in the background) as a function of TIME. Simulatenously, the gravitional three-body dynamics is evolved.
     
``` 
Evolve triple system

Evolve single stars (ignore entries for dummy secondary)


      TIME      M1       M2   K1 K2        SEP    ECC  R1/ROL1 R2/ROL2  TYPE
     0.0000   34.86861    0.00100  1 14 63792200.00000  0.00   0.000   0.000  INITIAL 
     5.8129   34.83353    0.00100  2 14 63856440.00000  0.00   0.000   0.000  KW_CHNGE
     5.8199   34.83309    0.00100  4 14 63857256.00000  0.00   0.000   0.000  KW_CHNGE
     6.3532   34.40681    0.00100  5 14 64648368.00000  0.00   0.000   0.000  KW_CHNGE
     6.3643   33.88068    0.00100 14 14 65666576.00000  0.01   0.000   0.000  KW_CHNGE
 15000.0000   33.88068    0.00100 14 14 65629124.00000  0.01   0.000   0.000  MAX_TIME
 Fallback:   1.0000000000000000     

      TIME      M1       M2   K1 K2        SEP    ECC  R1/ROL1 R2/ROL2  TYPE
     0.0000   31.64016    0.00100  1 14 61759344.00000  0.00   0.000   0.000  INITIAL 
     6.2848   31.61115    0.00100  2 14 61816020.00000  0.00   0.000   0.000  KW_CHNGE
     6.2930   31.61073    0.00100  4 14 61816848.00000  0.00   0.000   0.000  KW_CHNGE
     6.8721   31.33694    0.00100  5 14 62356908.00000  0.00   0.000   0.000  KW_CHNGE
     6.8841   29.58371    0.00100 14 14   155942.20312 99.99   0.000  -2.000  DISRUPT 
 15000.0000   29.58371    0.00100 14 14        0.00000 -1.00   0.000   0.000  MAX_TIME
 Fallback:  0.95997883689236163     

      TIME      M1       M2   K1 K2        SEP    ECC  R1/ROL1 R2/ROL2  TYPE
     0.0000   24.50586    0.00100  1 14 56717156.00000  0.00   0.000   0.000  INITIAL 
     7.9471   24.48971    0.00100  2 14 56754564.00000  0.00   0.000   0.000  KW_CHNGE
     7.9600   24.48932    0.00100  4 14 56755456.00000  0.00   0.000   0.000  KW_CHNGE
     8.7027   24.29811    0.00100  5 14 57202060.00000  0.00   0.000   0.000  KW_CHNGE
     8.7178   14.50907    0.00100 14 14      691.86182 99.99   0.000  -2.000  DISRUPT 
 15000.0000   14.50907    0.00100 14 14        0.00000 -1.00   0.000   0.000  MAX_TIME
 Fallback:  0.59707766402046147
```
After around 6.4 Myr an event was detected which deserves further attention. Here, the primary blew off in a supernova. The code updates the orbital parameters by taking into account any supernova kicks. In that case, the system remained stable and the integration is continued ... 

```
A termination event occured.

Primary supernova at 6.3643446
Plotting...
Plot saved as ./../plots/139_00001.png
Both orbits remain bound.
```

... but the supernova of the secondary disrupted the system and the program ends.

```
A termination event occured.

Secondary supernova at 6.88413429
Plotting...
Plot saved as ./../plots/139_00001.png
Inner orbit gets unbound.

```

The total evolution of the system is plotted in the directory `plots` and stored as a csv table in the directory `data`.

### Parallel runs

If you consider running a large population in parallel, consider using `gnu parallel`, e.g.,

`parallel -j 12 "python tse.py --random TRUE --seed {1} --Z {2}" ::: {1..1000} ::: 0.0002 0.002 0.02`

which would employ `-j 12` cores to evolve 1000 systems at three different metallicities (0.0002, 0.002, 0.02) in parallel. Note that it is not supported to run other parameters than seed and metallicity in parallel. If you want to explore other parameters you need to run them one after the other, e.g.,

`parallel -j 12 "python tse.py --alpha1 1.0 --random TRUE --seed {1} --Z {2}" ::: {1..1000} ::: 0.0002 0.002 0.02`

and then 

`parallel -j 12 "python tse.py --alpha1 2.0 --random TRUE --seed {1} --Z {2}" ::: {1..1000} ::: 0.0002 0.002 0.02`

and **not**

`parallel -j 12 "python tse.py --random TRUE --seed {1} --Z {2} --alpha1 {3}" ::: {1..1000} ::: 0.0002 0.002 0.02 ::: 1.0 2.0`. 

Otherwise, different threads will use the same mobse input and output files...

`src/serial_batch.sh` and `src/runtask` might be helpful scripts if you want to run the code on an HPC cluster that operates with slurm.


### Integration variables

y[0:3] : Eccentricity vector of the inner binary

y[3:6] : Dimensionless angular momentum vector of the inner binary

y[6] : Semi-major axis (Rsun) of the inner binary

y[7:10] : Eccentricity vector of the outer binary

y[10:13] : Dimensionless angular momentum vector of the outer binary

y[13] : Semi-major axis (Rsun) of the outer binary

y[14:17] : Rotation vector of the primary star (1/Myr), only evolved as long as the star is not a compact object

y[17:20] : Rotation vector of the secondary star (1/Myr), only evolved as long as the star is not a compact object

y[20:23] : BH spin vector of the primary, only evolved when the primary is a BH (norm 1, i.e., maximally spinning)

y[23:26] : BH spin vector of the secondary, only evolved when the secondary is a BH (norm 1, i.e., maximally spinning)

### Citation

When using this code for scientific publications cite

*Stegmann J., Antonini F. & Moe M., Mon.Not.Roy.Astron.Soc. 516 (2022) 1, 1406-1427* (https://doi.org/10.1093/mnras/stac2192)

along with the papers about mobse

*Giacobbo N. & Mapelli M., Mon.Not.Roy.Astron.Soc. (2018) 480, 2011–2030* (https://doi.org/10.1093/mnras/sty1999)

*Giacobbo N., Mapelli M., Spera M., Mon.Not.Roy.Astron.Soc. (2018) 474, 2959–2974* (https://doi.org/10.1093/mnras/stx2933)
