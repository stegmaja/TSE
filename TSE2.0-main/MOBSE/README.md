************************************************************************
```
    ╔╦╗╔═╗╔╗ ╔═╗╔═╗  
    ║║║║ ║╠╩╗╚═╗║╣   
    ╩ ╩╚═╝╚═╝╚═╝╚═╝  
```
(Massive Objects in Binary Stellar Evolution) 
===

#### Updated and customized version of BSE implemented by Nicola Giacobbo, in collaboration with Michela Mapelli and Mario Spera. 
**MOBSE** is a customized version of [BSE](http://astronomy.swin.edu.au/~jhurley/) (acronym for Binary Stellar Evolution], [Hurley et al. 2002](https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract)). With respect to BSE, MOBSE is a population-synthesis code that includes some important upgrades for the evolution of massive stars both singles and in binary systems.
Information on the updates present in **MOBSE** can be found at this [link](https://mobse-webpage.netlify.com/) and in the following papers: 

*"Revising Natal Kick Prescriptions in Population Synthesis Simulations"* [Giacobbo N. & Mapelli M., ApJ, Vol, 891, (2020)](https://ui.adsabs.harvard.edu/abs/2020ApJ...891..141G/abstract)  
*"The impact of electron-capture supernovae on merging double neutron stars"*
 [Giacobbo N. & Mapelli M., MNRAS 482, 2234–2243 (2019)](https://ui.adsabs.harvard.edu/abs/2019MNRAS.482.2234G/abstract)  
*"The progenitors of compact-object binaries: impact of metallicity, common envelope and natal kicks"*
 [Giacobbo N. & Mapelli M., MNRAS 480, 2011–2030 (2018)](https://ui.adsabs.harvard.edu/abs/2018MNRAS.480.2011G/abstract)  
*"Merging black hole binaries: the effects of progenitor's  metallicity, mass-loss rate and Eddington factor"*
 [Giacobbo N., Mapelli M., Spera M., MNRAS 474, 2959–2974 (2018)](https://ui.adsabs.harvard.edu/abs/2018MNRAS.474.2959G/abstract)  
*"The cosmic merger rate of stellar black hole binaries from the Illustris simulation"* [Mapelli et al., MNRAS 472, 2422-2435 (2017)](https://ui.adsabs.harvard.edu/abs/2017MNRAS.472.2422M/abstract) 

Any queries that are not answered by referring to these texts, or by 
reading the comments in the programs, can be addressed to: 
  *giacobbo.nicola@gmail.com*.
  
**In case you publish results based on this code, you should cite the following papers**:  
[Giacobbo N. & Mapelli M., MNRAS 480, 2011–2030 (2018)](https://ui.adsabs.harvard.edu/abs/2018MNRAS.480.2011G/abstract)  
[Giacobbo N., Mapelli M., Spera M., MNRAS 474, 2959–2974 (2018)](https://ui.adsabs.harvard.edu/abs/2018MNRAS.474.2959G/abstract)  

************************************************************************
Description of the MOBSE package
====

The MOBSE package contains the following FORTRAN files grouped in two folders (**src** and **input**): 

**IMORTANT**: many subroutines are named as in the BSE package but contain some 
          upgrades. They are identified by (#). There are some new 
          functions indicated by (+).  

In **src**:  
**mobse.f**       -> Main routine. Evolves one binary and creates data files (the
                    corrisponding bse.f)  
**evolve.f**      -> routine that controls the evolution of the binary (#/+ it 
                    is a modified version of evolv2.f)   
**comenv.f**      -> common envelope evolution (#)  
**corerd.f**      -> estimates the core radius of a giant-like star  
**deltat.f**      -> determines stellar evolution update timestep  
**dgcore.f**      -> determines the outcome when two degenerate cores merge  
**gntage.f**      -> calculates parameters of new star resulting from a merger  
**hrdiag.f**      -> decides which evolution stage the star is currently at 
                    and then the appropriate luminosity, radius
                    and core mass are calculated (#)   
**instar.f**      -> sets the collision matrix (#)  
**kick.f**        -> generates supernova kick and adjusts orbital parameters (#)  
**mix.f**         -> models stellar collisions (#)  
**mlwind.f**      -> contains the mass loss prescription (#)  
**mrenv.f**       -> calculates envelope parameters  
**ran3.f**        -> random number generator  
**rl.f**          -> calculates Roche-lobe radius   
**star.f**        -> derives the landmark timescales and luminosities 
                    that divide the various evolution stages (#)  
**zcnsts.f**      -> sets all the constants of the formulae which depend on 
                    metallicity (apart for stellar winds) (#)  
**zfuncs.f**      -> all the formulae as a collection of separate functions (#)  
**pisn.f**        -> contains the prescriptions for the pair-instability and the 
                    pulsation-pair-instability (+)  
**eddington.f**   -> calculates the Eddington factor (+)  
**fallback.f**    -> computes the fallback factor (+)  
**Makefile**      -> gfortran compiler. Use command `make mobse`

In **input**:  
**const_mobse.h** -> parameter file (#/+) and it contains new parameters  
**zdata.h**       -> contains all the coefficient values for zcnsts (#)  
**binary.input**  -> input file for mobse.f (#)  

and 

************************************************************************

MOBSE works as BSE and can be use in the same way. Then if you are 
familiar with BSE, the main routine mobse.f is an example 
to show how EVOLVE should be used. 
If have never used BSE, don't worry. The routines contain all original 
information and all the upgrades are commented. 

In the following I report a part of the README_BSE file which is valid also for
MOBSE, if you replace EVOLV2 with EVOLVE and bse.f with mobse.f

~~~~ 
In the case of EVOLV2 being called in a loop over many stars 
be careful to initialise each new star, i.e. 

mass(i) = mass0(i)
kstar(i) = 1
epoch(i) = 0.0
ospin(i) = 0.0
tphys = 0.0 

as well as setting the masses (mass0), period (tb) and eccentricity (ecc). 

However, the routine ZCNSTS only needs to be called each time you change
metallicity.

You may not want to use bse.f at all and instead, for example, prefer to use 
EVOLV2 directly as the main routine. 
Also, you may want to utilize the individual subroutines in different ways.
 ~~~~

************************************************************************
Definitions of the various variables and arrays
===
       MASS    Initial stellar mass in solar units  
       AJ      Current age in Myr  
       MT      Current mass in solar units (used for R)  
       TM      Main sequence lifetime  
       TN      Nuclear burning lifetime assuming no further mass loss  
       TSCLS   Evolution timescales for different stages  
       LUMS    Landmark luminosities  
       GB      Giant branch parameters  
       ZPARS   Parameters for distinguishing various mass intervals  
       R       Stellar radius in solar units  
       KW      Evolution type (0 - 15)  
       MC      Core mass in solar units  

### Evolution types for the stars (KW):  

              0 - deeply or fully convective low mass MS star  
              1 - Main Sequence star  
              2 - Hertzsprung Gap  
              3 - First Giant Branch  
              4 - Horizontal Branch / Core Helium Burning  
              5 - First Asymptotic Giant Branch / Red Supergiant  
              6 - Second Asymptotic Giant Branch  
              7 - Main Sequence Naked Helium star  
              8 - Hertzsprung Gap Naked Helium star  
              9 - Giant Branch Naked Helium star  
             10 - Helium White Dwarf  
             11 - Carbon/Oxygen White Dwarf  
             12 - Oxygen/Neon White Dwarf  
             13 - Neutron Star  
             14 - Black Hole  
             15 - Massless Supernova  
************************************************************************
### Timescales in TSCLS (they are all in Myr units) 
              1; BGB              2; He ignition   3; He burning  
              4; Giant t(inf1)    5; Giant t(inf2) 6; Giant t(Mx)  
              7; FAGB t(inf1)     8; FAGB t(inf2)  9; FAGB  t(Mx)  
             10; SAGB t(inf1)    11; SAGB t(inf2) 12; SAGB  t(Mx)  
             13; TP              14; t(Mcmax)  
************************************************************************
### Luminosities in LUMS (all luminosities are in solar units) 
              1; ZAMS             2; End MS        3; BGB  
              4; He ignition      5; He burning    6; L(Mx)  
              7; BAGB             8; TP  
************************************************************************
### GB = giant branch parameters 
              1; effective A(H)   2; A(H,He)       3; B  
              4; D                5; p             6; q  
              7; Mx               8; A(He)         9; Mc,BGB  
************************************************************************
### ZPARS = parameters that depends on the metallicity (Z)
              1; M below which hook doesn't appear on MS, Mhook.  
              2; M above which He ignition occurs non-degenerately, Mhef.  
              3; M above which He ignition occurs on the HG, Mfgb.  
              4; M below which C/O ignition doesn't occur, Mup.  
              5; M above which C ignites in the centre, Mec.  
              6; value of log D for M<= zpars(3)  
              7; value of x for Rgb propto M^(-x)  
              8; value of x for tMS = MAX(tHOOK,x*tBGB)  
              9; constant for McHeIf when computing Mc,BGB, mchefl.  
             10; constant for McHeIf when computing Mc,HeI, mchefl.  
             11; hydrogen abundance.  
             12; helium abundance.  
             13; constant x in rmin = rgb*x**y used by LM CHeB.  
             14; z**0.4 to be used for WD L formula.  
************************************************************************

Have fun!
