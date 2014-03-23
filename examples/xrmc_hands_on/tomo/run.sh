for th in $(seq 0 30 180); do
    cat quadric.templ | sed "s/-THETA-/$th/" > quadric.dat
    cat input.templ | sed "s/-THETA-/$th/" > input.dat
    xrmc input.dat
done
