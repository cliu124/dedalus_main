conda activate dedalus2
now="dedalus_$(date +"%Y%m%d%H%M%S")"
echo $now
mpiexec -n 1 python3 dedalus_main.py
mkdir -p "$now"
cp -r analysis "$now"
cp -r "$now" /mnt/d/Data/dedalus
echo $now
echo "finished"