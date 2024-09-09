F90 = gfortran	

CFLAG = -O2 -g

all: scvt_ini scvt_opt grid_ref draw_diag scvt.m comp_qual

scvt_ini: scvt_ini.f randgen.f density.f
	$(F90) $(CFLAG) -o scvt_ini scvt_ini.f randgen.f density.f 
scvt_opt: scvt_opt.f density.f process.f svtgen.f
	$(F90) $(CFLAG) -o scvt_opt scvt_opt.f density.f process.f \
		svtgen.f
grid_ref: grid_ref.f svtgen.f
	$(F90) $(CFLAG) -o grid_ref grid_ref.f svtgen.f
draw_diag: draw_diag.f svtgen.f
	$(F90) $(CFLAG) -o draw_diag draw_diag.f svtgen.f
comp_qual: comp_qual.f process.f svtgen.f density.f assist.f
	$(F90) $(CFLAG) -o comp_qual comp_qual.f process.f svtgen.f \
		density.f assist.f

clean:
	rm scvt_ini scvt_opt grid_ref draw_diag comp_qual

archive:
	tar -czvf SCVT.tgz scvt_ini.f scvt_opt.f grid_ref.f \
		draw_diag.f density.f process.f randgen.f assist.f\
		svtgen.f scvt.in scvt.m README Makefile comp_qual.f\
		scvt_uni_node_12.dat scvt_uni_node_42.dat \
		scvt_uni_node_162.dat scvt_uni_node_642.dat \
		scvt_uni_node_2562.dat scvt_uni_node_10242.dat
