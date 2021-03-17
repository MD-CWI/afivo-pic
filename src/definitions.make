OBJS := m_config.o m_geometry.o m_field.o m_domain.o m_refine.o m_time_step.o	\
m_particles.o m_globals.o m_user_methods.o m_photons.o m_init_cond.o

# Dependency information
m_domain.o: m_config.mod
m_domain.o: m_globals.mod
m_field.o: m_domain.mod
m_field.o: m_geometry.mod
m_field.o: m_globals.mod
m_field.o: m_user_methods.mod
m_globals.o: m_config.mod
m_init_cond.o: m_config.mod
m_init_cond.o: m_geometry.mod
m_init_cond.o: m_globals.mod
m_init_cond.o: m_user_methods.mod
m_particles.o: m_config.mod
m_particles.o: m_domain.mod
m_particles.o: m_field.mod
m_particles.o: m_globals.mod
m_particles.o: m_photons.mod
m_photons.o: m_config.mod
m_photons.o: m_domain.mod
m_photons.o: m_globals.mod
m_photons.o: m_time_step.mod
m_refine.o: m_config.mod
m_refine.o: m_geometry.mod
m_refine.o: m_globals.mod
m_time_step.o: m_config.mod
m_time_step.o: m_globals.mod
