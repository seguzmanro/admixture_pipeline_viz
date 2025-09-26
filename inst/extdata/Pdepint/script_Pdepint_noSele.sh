for k in $(seq 2 13)
do
admixture Pdepint_M095_noSele.bed $k --cv -j12 -B1000 | tee Pdepint_log${k}.out
done
