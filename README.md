# TSE

## Triple stellar evolution code

Single stellar evolution and binary interactions are modelled using the binary population synthesis code MOBSE (https://mobse-webpage.netlify.app) which builds upon the code BSE (https://www.ascl.net/1303.014). A copy of MOBSE is added to this repository. We have slighlty modified the files `mobse/src/evolve.f` and `mobse/src/mobse.f`. These changes were only made to ensure compatibility between TSE and MOBSE, e.g., change of array lengths and MOBSE input/output, but the physical prescription of MOBSE remains unaltered.  Using other versions of MOBSE or BSE-type codes is in general possible, as long as both modified files are maintained.

TSE evolves the secular dynamics of hierarchical triples and combines it with the stellar evolution of the inner binary and the tertiary companion.

### Installation

In order to install the code follow these steps to install

1. Navigate to the `mobse/src` directory.
2. Execute the `Makefile` by running `make mobse`.
3. Navigate back to the root directory.
4. Create two directories where the outcome of the simulation is stored, `mkdir plots` and `mkdir data`.
5. Create a conda environment with all the required packages by executing `conda env create -f environment.yml`.
6. Activate the newly created environment by `conda activate tse`.

Alternatively, simply execute `compile.sh` in the root directory. You may want to change the `mobse/src/Makefile` according to the specifications of your fortran compiler.

### Running

To run a system activate the environment (see above) and go to the `src` directory. Then execute

`python main.py`.

You may want to check out the number of optional arguments that are available by running

`python main.py -h`.

E.g., do 

`python main.py --random TRUE --seed 1 --Z 0.02`

to evolve a triple from a random distribution at solar metallicity.

If you are looking whether a triple leads to a particular configuration (e.g., the formation of BH triples) check out the function `CustomEvent` in `src/main.py` (including a working example).

TSE runs stably with the python version and packages listed in `requirements.txt`.

### Parallel runs

If you consider running a large population in parallel, consider using `gnu parallel`, e.g.,

`parallel -j 12 "python main.py --random TRUE --seed {1} --Z {2}" ::: {1..1000} ::: 0.0002 0.002 0.02`

which would employ `-j 12` cores to evolve 1000 systems at three different metallicities (0.0002, 0.002, 0.02) in parallel. Note that it is not supported to run other parameters than seed and metallicity in parallel. If you want to explore other parameters you need to run them one after the other, e.g.,

`parallel -j 12 "python main.py --alpha1 1.0 --random TRUE --seed {1} --Z {2}" ::: {1..1000} ::: 0.0002 0.002 0.02`

and then 

`parallel -j 12 "python main.py --alpha1 2.0 --random TRUE --seed {1} --Z {2}" ::: {1..1000} ::: 0.0002 0.002 0.02`

and **not**

`parallel -j 12 "python main.py --random TRUE --seed {1} --Z {2} --alpha1 {3}" ::: {1..1000} ::: 0.0002 0.002 0.02 ::: 1.0 2.0`. 

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
