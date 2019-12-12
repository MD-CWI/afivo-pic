ifndef NDIM
$(error NDIM is not set)
endif

ifndef MAIN_DIR
$(error MAIN_DIR is not set)
endif

AFIVO_DIR := $(MAIN_DIR)/afivo
PCORE_DIR := $(MAIN_DIR)/particle_core
LIBDIRS := $(MAIN_DIR)/lib_$(NDIM)d $(AFIVO_DIR)/lib_$(NDIM)d	\
$(AFIVO_DIR)/external_libraries/silo/lib $(AFIVO_DIR)/external_libraries/hypre/lib $(PCORE_DIR)
INCDIRS := $(TARGET_DIR) $(MAIN_DIR)/lib_$(NDIM)d $(AFIVO_DIR)/lib_$(NDIM)d $(PCORE_DIR)
LIBS := apic afivo silo HYPRE particle_core
PROG := apic

.PHONY: all clean allclean always_recompile

all: $(PROG)

clean:
	$(RM) *.o *.mod $(PROG)

allclean: clean
	$(MAKE) -C $(MAIN_DIR)/lib_$(NDIM)d clean

vpath %.f90 $(MAIN_DIR)/src

# Include compilation rules
include  $(MAIN_DIR)/makefiles/makerules.make

# Optionally include a local makefile
-include local.make

FFLAGS += -DNDIM=$(NDIM)

$(PROG): $(PROG).f90
	$(FC) -o $@ $(filter %.f90 %.o, $^) $(FFLAGS) $(addprefix -I,$(INCDIRS)) \
	$(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

# Dependencies
$(PROG): $(MAIN_DIR)/lib_$(NDIM)d/libapic.a
$(PROG): m_user.o

m_user.o: $(MAIN_DIR)/lib_$(NDIM)d/libapic.a

$(MAIN_DIR)/lib_$(NDIM)d/libapic.a: always_recompile
	$(MAKE) -C $(MAIN_DIR)/lib_$(NDIM)d

