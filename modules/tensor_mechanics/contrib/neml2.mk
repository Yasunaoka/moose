# There are 2 influential conf vars.
# They are controlled by the MOOSE configure script.
#
# ENABLE_LIBTORCH: True if we should compile with libTorch AND an installation
#                  of libTorch is available. False otherwise.
# LIBTORCH_DIR:    The path to the libTorch installation. NEML2 will use this
#                  libTorch as the backend.
#
# There are 2 influential env vars.
#
# ENABLE_NEML2:  If true, MOOSE will attempt to compile with NEML2 at
#                the NEML2_DIR. If false, MOOSE will not be compiled with
#                NEML2 (even if the NEML2 submodule is init'ed or a valid
#                NEML2_DIR is supplied). This variable defaults to true.
# NEML2_DIR:     The path to a valid NEML2 checkout. This variable defaults to
#                the NEML2 submodule. The user is responsible for setting this
#                variable if a custom NEML2 should be used. The compilation
#                will terminate with an error message if ENABLE_NEML2 is true
#                AND a NEML2 checkout cannot be found.

NEML2_DIR ?= $(MOOSE_DIR)/modules/tensor_mechanics/contrib/neml2
ifeq ($(wildcard $(NEML2_DIR)/CMakeLists.txt),)
_ENABLE_NEML2_DEFAULT = false
else
ifneq ($(ENABLE_LIBTORCH),true)
_ENABLE_NEML2_DEFAULT = false
else
_ENABLE_NEML2_DEFAULT = true
endif
endif

ifdef ENABLE_NEML2
_ENABLE_NEML2_SET_BY_USER = true
else
ENABLE_NEML2 = $(_ENABLE_NEML2_DEFAULT)
_ENABLE_NEML2_SET_BY_USER = false
endif

ifneq ($(ENABLE_NEML2),true)

ifeq ($(_ENABLE_NEML2_SET_BY_USER),true)
$(info Not compiling MOOSE with NEML2 because ENABLE_NEML2 is set to false.)
else
$(info Not compiling MOOSE with NEML2 because NEML2_DIR is not a valid NEML2 checkout and/or libtorch is not enabled.)
endif

else

ifneq ($(ENABLE_LIBTORCH),true)
$(error Attempting to compile MOOSE with NEML2, but libTorch is not enabled. \
  To enable libTorch, configure MOOSE with the --with-libtorch option. \
  To disable NEML2, set ENABLE_NEML2 to false)
endif

ifeq ($(wildcard $(NEML2_DIR)/CMakeLists.txt),)
ENABLE_NEML2 = false
ifeq ($(_ENABLE_NEML2_SET_BY_USER),true)
$(error Not compiling MOOSE with NEML2 because a valid NEML2 checkout cannot be found.  \
  To use the default NEML2 that comes with MOOSE, run `unset NEML2_DIR` and `git submodule update --init contrib/neml2`. \
	To use a custom NEML2, set the environment variable NEML2_DIR to an appropriate path. \
	To disable NEML2, set ENABLE_NEML2 to false.)
else
$(info Not compiling MOOSE with NEML2 because a valid NEML2 checkout cannot be found.  \
  To use the default NEML2 that comes with MOOSE, run `unset NEML2_DIR` and `git submodule update --init contrib/neml2`. \
	To use a custom NEML2, set the environment variable NEML2_DIR to an appropriate path. \
	To suppress this warning, set ENABLE_NEML2 to false.)
endif
endif

endif


###############################################################################
# At this point, we have everything needed to compile MOOSE with NEML2
###############################################################################
ifeq ($(ENABLE_NEML2),true)

$(info Compiling MOOSE with NEML2.)

app_non_unity_dirs += $(shell find $(APPLICATION_DIR)/src/neml2 -type d -not -path '*/.libs*' 2> /dev/null)
app_non_unity_dirs += $(shell find $(APPLICATION_DIR)/test/src/neml2 -type d -not -path '*/.libs*' 2> /dev/null)

NEML2_INCLUDE        := $(NEML2_DIR)/include
NEML2_SRC            := $(shell find $(NEML2_DIR)/src -name "*.cxx")
NEML2_OBJ            := $(patsubst %.cxx,%.$(obj-suffix),$(NEML2_SRC))
NEML2_LIB            := $(NEML2_DIR)/libNEML2-$(METHOD).la

$(NEML2_LIB): $(NEML2_OBJ)
	@echo "Linking Library "$@"..."
	@$(libmesh_LIBTOOL) --tag=CC $(LIBTOOLFLAGS) --mode=link --quiet \
	  $(libmesh_CC) $(libmesh_CFLAGS) -o $@ $(NEML2_OBJ) $(libmesh_LDFLAGS) $(EXTERNAL_FLAGS) -rpath $(NEML2_DIR)
	@$(libmesh_LIBTOOL) --mode=install --quiet install -c $(NEML2_LIB) $(NEML2_DIR)

$(NEML2_DIR)/src/%.$(obj-suffix) : $(NEML2_DIR)/src/%.cxx
	@echo "Compiling C++ (in "$(METHOD)" mode) "$<"..."
	@$(libmesh_LIBTOOL) --tag=CXX $(LIBTOOLFLAGS) --mode=compile --quiet \
	  $(libmesh_CXX) $(libmesh_CPPFLAGS) $(ADDITIONAL_CPPFLAGS) $(libmesh_CXXFLAGS) $(app_INCLUDES) $(libmesh_INCLUDE) -w -DHAVE_CONFIG_H -MMD -MP -MF $@.d -MT $@ -c $< -o $@

ADDITIONAL_INCLUDES  += -iquote$(NEML2_INCLUDE)
ADDITIONAL_CPPFLAGS  += -DNEML2_ENABLED -DDTYPE=Float64
ADDITIONAL_LIBS      += -L$(NEML2_DIR) -lNEML2-$(METHOD)
ADDITIONAL_DEPEND_LIBS += $(NEML2_LIB)

endif
