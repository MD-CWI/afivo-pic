INCDIRS	:= afivo/src particle_core/src
LIBDIRS := afivo/src afivo/external_libraries/silo/lib particle_core/src
LIBS	:= afivo silo particle_core
VPATH	:= config_fortran

include afivo/makerules.make

OBJS_XD := m_field_Xd.o m_init_cond_Xd.o m_photoi_Xd.o
OBJS := m_streamer.o m_config.o m_geometry.o m_transport_data.o
PROGS := ../apic_3d

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
../apic_2d: $(OBJS) $(OBJS_XD:%_Xd.o=%_2d.o)
../apic_3d: $(OBJS) $(OBJS_XD:%_Xd.o=%_3d.o)

apic_2d.o: m_streamer.mod $(OBJS_XD:%_Xd.o=%_2d.mod)
apic_3d.o: m_streamer.mod $(OBJS_XD:%_Xd.o=%_3d.mod)
m_streamer.o: m_config.mod m_transport_data.mod
m_init_cond_2d.o m_init_cond_3d.o: m_geometry.mod

