for irun in $(seq 1 10); do
    for ith in $(seq -90 10 90); do
	echo $irun " " $ith
	qsub -p 0.7 -V -q hpcserial job1.sh $PWD $irun $ith
	sleep 3
    done
done
