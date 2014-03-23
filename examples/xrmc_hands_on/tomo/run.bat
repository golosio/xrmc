for /l %i in (0, 10, 180); do cat quadric.templ | sed \"s/_THETA_/%i/\" > quadric.dat & cat input.templ | sed \"s/_THETA_/%i/\" > input.dat & xrmc input.dat 
