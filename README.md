# TSE

## Triple stellar evolution code

Single stellar evolution and binary interactions are modelled using the binary population synthesis code MOBSE (https://mobse-webpage.netlify.app) which builds upon the code BSE (https://www.ascl.net/1303.014). A copy of mobse is added to this repository. Changes are only made to the main file mobse.f to get all outputs which are used by TSE. Using other versions of MOBSE or BSE-type codes in general is possible, as long as the main file is replaced by this modified version of mobse.f.

TSE evolves the secular dynamics of hierarchical triples and combines it with the stellar evolution of the inner binary and the tertiary companion.

In order to run the code first compile mobse by running 

`make mobse`

in the mobse directory and

`mkdir plots`

in the root directory. Alternatively, simply execute `compile.sh` in the root directory. You may want to change the Makefile according to the specifications of your fortran compiler.

To run a system go to the src directory and run

`python main.py`.

You may want to check out the number of optional arguments that are available by running

`python main.py -h`.

E.g., do 

`python main.py --random TRUE --seed 1 --Z 0.02`

to evolve a triple from a random distribution at solar metallicity.

If you are looking whether a triples leads to a particular configuration (e.g., the formation of BH triples) check out the function CustomEvent in main.py (including a working example).
