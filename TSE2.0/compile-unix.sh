module load anaconda/2020.02/2019.03
rm *.so __pycache__/* *.mod
f2py -c MOBSE/input/const_mobse.h MOBSE/input/zdata.h -m my_lib MOBSE/src/comenv.f MOBSE/src/corerd.f MOBSE/src/deltat.f MOBSE/src/dgcore.f MOBSE/src/evolve.f MOBSE/src/gntage.f \
MOBSE/src/hrdiag.f MOBSE/src/instar.f MOBSE/src/kick.f MOBSE/src/mix.f MOBSE/src/mlwind.f MOBSE/src/mrenv.f MOBSE/src/ran3.f MOBSE/src/rl.f \
MOBSE/src/star.f MOBSE/src/zcnsts.f MOBSE/src/zfuncs.f MOBSE/src/pisn.f MOBSE/src/eddington.f MOBSE/src/fallback.f