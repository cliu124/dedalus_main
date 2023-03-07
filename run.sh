conda activate dedalus2
now=$(date +"%Y%m%d%H%M%S")
WORKDIR= dedalus_$now
mpiexec -n 1 python3 dedalus_main.py
mkdir "$WORKDIR"
cp -r analysis "$WORKDIR"
cp -r "$WORKDIR" /mnt/d/Data/dedalus
echo "$WORKDIR"
echo "finished"