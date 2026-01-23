# ====================================================================================
# Circular Tilting Tank Simulation Makefile
# ====================================================================================
#
# USAGE:
#   make [target] [OVERRIDE_VAR=value]
#
# TARGETS:
#   all       : Generate mesh, motion file, and setup case (default).
#   mesh      : Generate geometry (.geo) and mesh (.msh) using Gmsh.
#   motion    : Generate static 6DoF file + tilted gravity (constant/6DoF.dat, constant/g).
#   case      : Initialize field variables (runs update_setFields.py).
#   run       : Run the OpenFOAM simulation (gmshToFoam -> setFields -> foamRun).
#   clean     : Remove generated mesh, motion, and result files.
#   help      : Show this help message.
#
# PARAMETERS:
#   H             : Tank height (default: 1.0)
#   D             : Tank diameter (default: 0.5)
#   MESH_SIZE     : Mesh characteristic length (default: 0.05)
#   GEO_TYPE      : Geometry type: 'flat' or 'cap' (default: flat)
#   TILT_DEG      : Static tilt angle about +X axis in degrees (default: 5.0)
#   DURATION      : Simulation duration in seconds (default: 10.0)
#   DT            : Time step for motion file (default: 0.01)
#   RAMP_DURATION : Smooth start ramp duration (default: -1, i.e., 10% of total)
#
# ====================================================================================

# Parameters
H ?= 1.0
D ?= 0.5
MESH_SIZE ?= 0.05
TILT_DEG ?= 5.0
DURATION ?= 10.0
DT ?= 0.01
GEO_TYPE ?= flat

.PHONY: all mesh motion case run clean help

all: prep run

help:
	@grep -E '^[a-zA-Z_-]+:.*?$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-15s\033[0m %s\n", $$1, $$2}'
	@echo ""
	@echo "Read the header in Makefile for variable documentation."

# Prepare geometry and coordinate files (Run on Host/Login Node)
prep: motion case
	python3 generate_mesh.py $(H) $(D) $(MESH_SIZE) $(GEO_TYPE)
	gmsh -3 cylinder.geo -format msh2 -o cylinder.msh

# Generate Motion
motion:
	python3 generate_tilt.py $(TILT_DEG) $(DURATION) $(DT)

# Configure Case (Initial Fields)
case:
	# Update setFieldsDict for half-fill height (H/2)
	python3 update_setFields.py $(H)

# Run Simulation
ifeq ($(OSCAR),1)
    OF_PREFIX := of13 
endif

OF_PREFIX ?= 
N_CPUS ?= 1
ADAPTIVE_STOP ?= 1
run:
	# commands to be run in OpenFOAM environment
	$(OF_PREFIX) gmshToFoam cylinder.msh
	$(OF_PREFIX) setFields
	@if [ $(N_CPUS) -gt 1 ]; then \
		echo "Running in parallel on $(N_CPUS) cores..."; \
		$(OF_PREFIX) decomposePar -force; \
		if [ $(ADAPTIVE_STOP) -eq 1 ]; then \
			$(OF_PREFIX) python3 adaptive_stop.py --parallel --np $(N_CPUS); \
		else \
			$(OF_PREFIX) mpirun -np $(N_CPUS) foamRun -parallel; \
		fi; \
		echo "Reconstructing parallel results..."; \
		$(OF_PREFIX) reconstructPar; \
		echo "Cleaning up processor directories..."; \
		rm -rf processor*; \
	else \
		echo "Running in serial..."; \
		if [ $(ADAPTIVE_STOP) -eq 1 ]; then \
			$(OF_PREFIX) python3 adaptive_stop.py; \
		else \
			$(OF_PREFIX) foamRun; \
		fi; \
	fi

resume:
	@if [ $(N_CPUS) -gt 1 ]; then \
		echo "Resuming in parallel on $(N_CPUS) cores..."; \
		if [ $(ADAPTIVE_STOP) -eq 1 ]; then \
			$(OF_PREFIX) python3 adaptive_stop.py --parallel --np $(N_CPUS); \
		else \
			$(OF_PREFIX) mpirun -np $(N_CPUS) foamRun -parallel; \
		fi; \
		echo "Reconstructing parallel results..."; \
		$(OF_PREFIX) reconstructPar; \
		echo "Cleaning up processor directories..."; \
		rm -rf processor*; \
	else \
		echo "Resuming in serial..."; \
		if [ $(ADAPTIVE_STOP) -eq 1 ]; then \
			$(OF_PREFIX) python3 adaptive_stop.py; \
		else \
			$(OF_PREFIX) foamRun; \
		fi; \
	fi

clean:
	rm -f cylinder.geo cylinder.msh
	rm -f constant/6DoF.dat
	rm -rf 0/alpha.water.orig  # created by setFields usually backup
	# standard cleanup
