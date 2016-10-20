N=10
python cohesion_random_generator.py $N > cohesion.dat
python phi_random_generator.py $N > phi.dat
COHESIONLIST=`cat cohesion.dat`
PHILIST=`cat phi.dat`
for C in $COHESIONLIST
do
	for PHI in $PHILIST
	do
		echo "Computing for COHESION=$C, PHI=$PHI,"
		./demo_wavetrial_pulseVM $C $PHI
	done	
done
