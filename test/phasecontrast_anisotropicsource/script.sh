for iE in $(seq 0 200); do
    E=$(xraylib "10. + $iE/10")
    A=$(xraylib "exp(-($iE/10 - 10)*($iE/10 - 10)/(2.*4*4))")
    echo $E $A
done
