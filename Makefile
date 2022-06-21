#------------------------------------
# Makefile 3D reconnection code on gateway
#------------------------------------
MPICH=

#ifeq ($(COMP),xl)
#F90= mpixlf
#EXE= ../bin/xl_turbo_rec.x
#else
F90= mpipgifort
EXE= ./pgi_turbo_rec.x
#endif

FFLAGS =  -O2 -g -c -acc=gpu -ta=tesla -Minfo=accel -Mcuda
LNKFLAGS = -acc=gpu -ta=tesla -Minfo=accel -I/cineca/prod/opt/compilers/hpc-sdk/2021/binary/Linux_ppc64le/21.5/cuda/11.3/include -L/cineca/prod/opt/compilers/hpc-sdk/2021/binary/Linux_ppc64le/21.5/cuda/11.3/lib64 -cudalib=cusparse -lnvToolsExt

all: $(COMM) MAIN

# definit une macro OBJ pour les fichiers objets

OBJ= nvtx.o gpu_func.o corrente_res.o der1x_free.o der1x_y.o \
     der1y.o der1z.o der2y.o filterx_t.o filtery_t.o filterz_t.o \
     filter_yz_mod.o ftrout.o gpmio.o helm_mod.o init_mod.o new_recon.o out_field.o \
     out_field_rhs.o out_grid.o out_inv.o outwrite.o poisson_mod.o pstartup.o \
     rhs_xy.o rhs_z.o trasponi_yx.o trasponi_zx.o tridLU.o \
     p01abft.o p01abzt.o s11aaft.o s11acft.o condinit_kh.o \
     x04aaft.o x04baft.o energie_mod.o
#heat_flux_mod.o condinit_mod_ran.o  x04aaft.o x04baft.o
#condinit_mod_ran_2D.o
#condinit_mod_ran_w3D.o
#condinit_kh.o

COMM= comm/comm.a

$(COMM):
	( cd comm ; $(F90) $(FFLAGS) *.f )
	( cd comm ; ar  rv comm.a *.o )


#$(FFT):
#       ( cd utility ; $(F90) $(FFLAGS) tempert.f )

# 1ere dependence: l'ex. d'un fichier .o depend de l'ex. d'un fichier .f
.f.o:
	$(F90) $(FFLAGS) $*.f

# 2eme dependence: Main, qui depend de la macro OBJ
MAIN:   $(OBJ) $(COMM) $(FFT)
	$(F90)  -o $(EXE) $(LNKFLAGS) $(OBJ) $(MODS) $(FFT) $(COMM) $(LIB)

clean:
	rm -f *.o *.mod *.a */*.o */*.mod */*.a *core */*core
