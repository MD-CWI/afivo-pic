# This is a template to build afivo-pic libraries in 2D/3D

ifndef NDIM
$(error NDIM is not set)
endif

ifndef MAIN_DIR
$(error MAIN_DIR is not set)
endif

AFIVO_DIR := $(MAIN_DIR)/afivo
AFIVO_LIB := $(AFIVO_DIR)/lib_$(NDIM)d/libafivo.a
PCORE_DIR := $(MAIN_DIR)/particle_core
PCORE_LIB := $(PCORE_DIR)/libparticle_core.a
LIB := libapic.a

 # Used in the compilation rules
INCDIRS := $(AFIVO_DIR)/lib_$(NDIM)d $(PCORE_DIR)

.PHONY: all clean always_recompile
all: $(LIB)

# Where to find the source files
vpath %.f90 $(MAIN_DIR)/src $(MAIN_DIR)/src/config_fortran \
	$(MAIN_DIR)/src/lookup_table_fortran $(MAIN_DIR)/src/rng_fortran

# Dependencies
include $(MAIN_DIR)/src/definitions.make

# Compilation rules
include  $(MAIN_DIR)/makefiles/makerules.make

FFLAGS += -DNDIM=$(NDIM)

$(OBJS): $(AFIVO_LIB) $(PCORE_LIB)

$(LIB): $(OBJS)
	$(RM) $@
	$(AR) rcs $@ $^

clean:
	$(RM) *.o *.mod $(LIB)

$(LIB): $(AFIVO_LIB) $(PCORE_LIB)

$(AFIVO_LIB): always_recompile
	$(MAKE) -C $(AFIVO_DIR) lib_$(NDIM)d

$(PCORE_LIB): always_recompile
	$(MAKE) -C $(PCORE_DIR)
