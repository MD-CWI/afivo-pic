INCDIRS	:= afivo/src particle_core/src
LIBDIRS := afivo/src afivo/external_libraries/silo/lib particle_core/src
LIBS	:= afivo silo particle_core
VPATH	:= config_fortran

include afivo/makerules.make

OBJS := apic_2d.o m_streamer.o m_field_2d.o m_config.o m_init_cond_2d.o	\
m_geometry.o m_transport_data.o m_photoi_2d.o
PROGS := ../apic_2d

.PHONY: all clean afivo/src/libafivo.a particle_core/src/libparticle_core.a

all: 	$(PROGS)

libstreamer.a: $(OBJS)
	$(RM) $@
	$(AR) rcs $@ $^

clean:
	$(RM) $(PROGS) *.o *.mod libstreamer.a

%_2d.f90: %_Xd.f90
	@echo "*** $< has changed, regenerating $@ ***"
	@$(PPROCESS_2D)

%_3d.f90: %_Xd.f90
	@echo "*** $< has changed, regenerating $@ ***"
	@$(PPROCESS_3D)

# Keep intermediate files for inspection
.PRECIOUS: %_2d.f90 %_3d.f90

# How to get executables from .o object files
../%: %.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

afivo/src/libafivo.a:
	$(MAKE) -C afivo

particle_core/src/libparticle_core.a:
	$(MAKE) -C particle_core

# Dependency information
$(PROGS): afivo/src/libafivo.a particle_core/src/libparticle_core.a
$(PROGS): $(OBJS)

apic_2d.o: m_streamer.mod m_field_2d.mod m_init_cond_2d.mod m_photoi_2d.mod
m_streamer.o: m_config.mod m_transport_data.mod
m_init_cond_2d.o: m_geometry.mod
