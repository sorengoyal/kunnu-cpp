N=10
python rho_random_generator.py $N > rho.dat
python e_random_generator.py $N > e.dat
RHOLIST=`cat rho.dat`
ELIST=`cat e.dat`
for E in $ELIST
do
	for RHO in $RHOLIST
	do
		echo "Computing for RHO=$RHO, E=$E"
		python 1material_wave.py $RHO $E
	done	
done
