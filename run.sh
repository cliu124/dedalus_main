conda activate dedalus2
mpiexec -n 1 python3 dedalus_main.py
cp -r analysis /mnt/d/Data/dedalus/local
