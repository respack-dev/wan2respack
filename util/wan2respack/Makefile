FC = ifort

FFLAGS = -O2 -qopenmp -xHost -g -traceback -warn all -check all
#LIBBLAS = -lmkl_intel_lp64 -Wl,--start-group -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread -lm 
LIBBLAS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm 

#--------------------------------------------------------------------
.SUFFIXES:
.SUFFIXES: .o .F90 .F .f90 .c
# default suffix rules. Do not edit this area.
.F.o:
	${FC} ${F90FLAGS} ${FFLAGS} -c $<
.f90.o:
	${FC} ${F90FLAGS} ${FFLAGS} -c $<
.F90.o:
	${FC} ${F90FLAGS} ${FFLAGS} -c $<
#--------------------------------------------------------------------

TARGET = gen_mk.x gen_wan.x
OBJECTS = generate_sample_mk.o m_rd_dat_wfn.o util.o generate_wan.o

all : $(TARGET)

gen_mk.x : generate_sample_mk.o m_rd_dat_wfn.o util.o
	${FC} ${LDFLAGS} -o $@ generate_sample_mk.o m_rd_dat_wfn.o util.o ${LIBBLAS} 

gen_wan.x : generate_wan.o util.o
	${FC} ${LDFLAGS} -o $@ generate_wan.o util.o ${LIBBLAS} 

# gen_mk.x
generate_sample_mk.o: generate_sample_mk.f90 m_rd_dat_wfn.o util.o 
m_rd_dat_wfn.o: m_rd_dat_wfn.f90 util.o 

# gen_wan.x
generate_wan.o: generate_wan.f90 util.o

clean:
	${RM} ${TARGET} ${OBJECTS} *.mod
