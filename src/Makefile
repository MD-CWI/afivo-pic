INCDIRS	:= afivo/src particle_core/src
LIBDIRS := afivo/src afivo/external_libraries/silo/lib particle_core/src
LIBS	:= afivo silo particle_core
VPATH	:= config_fortran

include afivo/makerules.make

OBJS_XD := m_field_Xd.o m_init_cond_Xd.o m_photoi_Xd.o m_domain_Xd.o	\
m_refine_Xd.o m_time_step_Xd.o m_particles_Xd.o m_globals_Xd.o
OBJS := m_config.o m_geometry.o m_transport_data.o
PROGS := ../apic_2d ../apic_3d

.PHONY: all libs clean allclean

all: 	$(PROGS)

libstreamer.a: $(OBJS)
	$(RM) $@
	$(AR) rcs $@ $^

clean:
	$(RM) $(PROGS) *.o *.mod libstreamer.a

allclean: clean
	$(MAKE) -C afivo clean
	$(MAKE) -C particle_core clean

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

$(OBJS): afivo/src/libafivo.a particle_core/src/libparticle_core.a
$(OBJS_XD:%_Xd.o=%_2d.o): afivo/src/libafivo.a particle_core/src/libparticle_core.a
$(OBJS_XD:%_Xd.o=%_3d.o): afivo/src/libafivo.a particle_core/src/libparticle_core.a

# Put content of list_module_deps.sh below
apic_2d.o: m_domain_2d.mod
apic_2d.o: m_field_2d.mod
apic_2d.o: m_globals_2d.mod
apic_2d.o: m_init_cond_2d.mod
apic_2d.o: m_particles_2d.mod
apic_2d.o: m_photoi_2d.mod
apic_2d.o: m_refine_2d.mod
apic_2d.o: m_time_step_2d.mod
apic_3d.o: m_domain_3d.mod
apic_3d.o: m_field_3d.mod
apic_3d.o: m_globals_3d.mod
apic_3d.o: m_init_cond_3d.mod
apic_3d.o: m_particles_3d.mod
apic_3d.o: m_photoi_3d.mod
apic_3d.o: m_refine_3d.mod
apic_3d.o: m_time_step_3d.mod
m_domain_2d.o: m_config.mod
m_domain_3d.o: m_config.mod
m_field_2d.o: m_domain_2d.mod
m_field_2d.o: m_globals_2d.mod
m_field_3d.o: m_domain_3d.mod
m_field_3d.o: m_globals_3d.mod
m_globals_2d.o: m_config.mod
m_globals_3d.o: m_config.mod
m_init_cond_2d.o: m_domain_2d.mod
m_init_cond_2d.o: m_geometry.mod
m_init_cond_2d.o: m_globals_2d.mod
m_init_cond_3d.o: m_domain_3d.mod
m_init_cond_3d.o: m_geometry.mod
m_init_cond_3d.o: m_globals_3d.mod
m_particles_2d.o: m_config.mod
m_particles_2d.o: m_domain_2d.mod
m_particles_2d.o: m_globals_2d.mod
m_particles_2d.o: m_photoi_2d.mod
m_particles_3d.o: m_config.mod
m_particles_3d.o: m_domain_3d.mod
m_particles_3d.o: m_globals_3d.mod
m_particles_3d.o: m_photoi_3d.mod
m_photoi_2d.o: m_config.mod
m_photoi_2d.o: m_domain_2d.mod
m_photoi_2d.o: m_globals_2d.mod
m_photoi_3d.o: m_config.mod
m_photoi_3d.o: m_domain_3d.mod
m_photoi_3d.o: m_globals_3d.mod
m_refine_2d.o: m_config.mod
m_refine_2d.o: m_geometry.mod
m_refine_2d.o: m_globals_2d.mod
m_refine_2d.o: m_init_cond_2d.mod
m_refine_2d.o: m_transport_data.mod
m_refine_3d.o: m_config.mod
m_refine_3d.o: m_geometry.mod
m_refine_3d.o: m_globals_3d.mod
m_refine_3d.o: m_init_cond_3d.mod
m_refine_3d.o: m_transport_data.mod
m_time_step_2d.o: m_config.mod
m_time_step_2d.o: m_globals_2d.mod
m_time_step_3d.o: m_config.mod
m_time_step_3d.o: m_globals_3d.mod
