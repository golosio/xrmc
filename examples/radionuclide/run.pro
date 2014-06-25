a=xrmc_read_angio_image2('output.dat')
help,a
;b=read_spe('stainlesssteel_4.spe')
my_plot = plot(dindgen(2000)/50.0, total(a.data[1:4, * , 0 , 0],1) > 1,'r',/ylog)
;my_plot = plot(dindgen(2000)/40.0, b.data > 1,'b',/ylog,/overplot)



END
