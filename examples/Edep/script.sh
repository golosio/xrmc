mu_ph=$(xraylib '1.06*(0.0805*CS_Photo(1,50)+0.5999*CS_Photo(6,50)+0.3196*CS_Photo(8,50))')

mu_C=$(xraylib '1.06*(0.0805*CS_Compt(1,50)+0.5999*CS_Compt(6,50)+0.3196*CS_Compt(8,50))')

Nd=21600

mu_tot=$(xraylib '1.06*(0.0805*CS_Total(1,50)+0.5999*CS_Total(6,50)+0.3196*CS_Total(8,50))')

:>profile_ph_th.txt
for i in $(seq 0 399); do
    if [ $i -lt 34 ] || [ $i -gt 367 ]; then
	echo "$i 0" >> profile_ph_th.txt
    else
	r=$(echo "94.0+0.03*$i" | bc -l)
	N0=$(echo "$Nd*115.0*115.0/($r*$r)" | bc -l)
	#echo $N0
	d=$(echo "0.03*($i-34)" | bc -l)
	Psurv=$(xraylib "exp(-$mu_tot*$d)")
	#echo "$d $mu_tot $Psurv"
	N1=$(echo "$N0*$Psurv" | bc -l)
	Eph=$(echo "50.0*$N1*0.03*$mu_ph" | bc -l)
	echo "$i $Eph" >> profile_ph_th.txt
    fi
done

:>profile_C_th.txt
for i in $(seq 0 399); do
    if [ $i -lt 34 ] || [ $i -gt 367 ]; then
	echo "$i 0" >> profile_C_th.txt
    else
	r=$(echo "94.0+0.03*$i" | bc -l)
	N0=$(echo "$Nd*115.0*115.0/($r*$r)" | bc -l)
	#echo $N0
	d=$(echo "0.03*($i-34)" | bc -l)
	Psurv=$(xraylib "exp(-$mu_tot*$d)")
	#echo "$d $mu_tot $Psurv"
	N1=$(echo "$N0*$Psurv" | bc -l)
	EC=$(echo "4.151*$N1*0.03*$mu_C" | bc -l)
	echo "$i $EC" >> profile_C_th.txt
    fi
done


