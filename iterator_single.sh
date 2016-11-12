N=10
python cohesion_phi_random_generator.py $N > cohesion_phi.dat
LIST=`cat cohesion_phi.dat`
for L in $LIST
do
	echo "Computing for COHESION and PHI =$L,"
        echo "mpirun -4 demo_wave2DPulseVM $L"
	mpirun -n 4 demo_wavetrial_pulseVM $L
done
