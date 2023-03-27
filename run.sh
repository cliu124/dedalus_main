conda activate dedalus2
now="dedalus_$(date +"%Y%m%d%H%M%S")"
ln -s /home/cliu124/dedalus_main/dedalus_20230326162736/analysis/analysis_s1.h5 restart.h5
mpiexec -n 1 python3 dedalus_main.py
echo $now
mkdir -p "$now"
cp -r analysis* "$now"
cp -r "$now" /mnt/d/Data/dedalus
echo $now
echo "finished"
rm -rf restart.h5