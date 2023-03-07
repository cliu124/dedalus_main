conda activate dedalus2
now=$(date+"%Y%m%d%H%M%S")
WORKDIR= dedalus_$now
mkdir -p "$WORKDIR" && cp -r dedalus*.py "$WORKDIR" && cp -r dedalus.cfg "$WORKDIR" && cp submit_dedalus "$WORKDIR" && cd "$WORKDIR" || exit -1
mpiexec -n 1 python3 dedalus_main.py
cd ..
cp -r "$WORKDIR" /mnt/d/Data/dedalus
