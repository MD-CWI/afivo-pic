OBJS := m_config.o m_geometry.o m_transport_data.o m_field.o m_init_cond.o	\
m_photoi.o m_domain.o m_refine.o m_time_step.o m_particles.o m_globals.o

# Dependency information
apic.o: m_domain.mod
apic.o: m_field.mod
apic.o: m_globals.mod
apic.o: m_init_cond.mod
apic.o: m_particles.mod
apic.o: m_photoi.mod
apic.o: m_refine.mod
apic.o: m_time_step.mod
m_domain.o: m_config.mod
m_field.o: m_domain.mod
m_field.o: m_globals.mod
m_globals.o: m_config.mod
m_init_cond.o: m_domain.mod
m_init_cond.o: m_geometry.mod
m_init_cond.o: m_globals.mod
m_particles.o: m_config.mod
m_particles.o: m_domain.mod
m_particles.o: m_field.mod
m_particles.o: m_globals.mod
m_particles.o: m_photoi.mod
m_photoi.o: m_config.mod
m_photoi.o: m_domain.mod
m_photoi.o: m_globals.mod
m_refine.o: m_config.mod
m_refine.o: m_geometry.mod
m_refine.o: m_globals.mod
m_refine.o: m_init_cond.mod
m_refine.o: m_transport_data.mod
m_time_step.o: m_config.mod
m_time_step.o: m_globals.mod
