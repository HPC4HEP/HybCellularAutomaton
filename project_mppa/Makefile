# Needed for Ubuntu/Debian since sh points to dash instead of bash...
SHELL := /bin/bash

# board can be developer, tc2 or emb01
board ?= developer

K1_TOOLCHAIN_DIR := $(shell echo $${K1_TOOLCHAIN_DIR:-/usr/local/k1tools/})

SRC_DIR   := src
HOST_SRC_DIR := ${SRC_DIR}/host
K1_IO_SRC_DIR := ${SRC_DIR}/io
K1_CLUSTER_SRC_DIR := ${SRC_DIR}/cluster
INCLUDE_SRC_DIR := ${SRC_DIR}/include

K1_IO_APPLI_NAME := project_io
K1_CLUSTER_APPLI_NAME := project_cluster
HOST_APPLI_NAME := project_host

# emb01 is always a remote target
# tc2 and developer can be both local or remote (local by default)
ifeq (${board}, emb01)
remote := true
else
remote := false
endif

# Define weither we are using OSes or not. can be os/rtems/nodeos/bare
# If not present, default to os
system-name := os

# Common cflags and lflags (for all targets)
# These flags will be propagated through the hierarchy
cflags = -g -O2 -I${INCLUDE_SRC_DIR} 
lflags := -lm

# The mppa-bin target is used to generate multibinary for specified target
mppa-bin := multibin
multibin-name := project_mppa.mpk

# Define which binaries will composed multibin1
# Note that the first binary in the list defines the boot binary
multibin-objs := ${K1_IO_APPLI_NAME} ${K1_CLUSTER_APPLI_NAME}
#multibin-boot := ${K1_IO_APPLI_NAME}
# You can pass additionnal flags to "k1-create-multibinary" using the
# following command.
# for example, '-w ".*"' will active all tracepoints
ifneq (${board}, emb01)
multibin-flags := -w ".*"
endif

# Common flags for k1 binaries
k1-cflags := -mboard=$(board)

# Cluster binaries
cluster-bin := ${K1_CLUSTER_APPLI_NAME}

# Each *-srcs variable describe sources for the binaries
${K1_CLUSTER_APPLI_NAME}-srcs := ${K1_CLUSTER_SRC_DIR}/cluster_main.c
${K1_CLUSTER_APPLI_NAME}-cflags := -fopenmp
${K1_CLUSTER_APPLI_NAME}-lflags := -fopenmp -lmppaipc 

# The io-bin var is used to build groups of var
io-bin := ${K1_IO_APPLI_NAME}
# Common cflags for io sources
io-cflags := -fopenmp

${K1_IO_APPLI_NAME}-srcs := ${K1_IO_SRC_DIR}/io_main.c
# Flags can be specified per binary
${K1_IO_APPLI_NAME}-cflags := -fopenmp 
${K1_IO_APPLI_NAME}-lflags := -lmppaipc

# Same rule as the other but this time for host
# It define the host programs we will compile
host-bin := project_host
# Common flags for host programs
host-cflags :=
host-lflags := 

${HOST_APPLI_NAME}-srcs := ${HOST_SRC_DIR}/host_main.c
${HOST_APPLI_NAME}-cflags :=
${HOST_APPLI_NAME}-lflags :=

# Hooks specified in help-hooks are rules which will be called after displaying
# the default help
help-hooks := my_help_rule
# Post build hooks will be called after build
post-build-hooks := post_build_rule
# clean hooks will be called with clean
clean-hooks := clean-traces

# Finally, include the Kalray Makefile which will instantiate the rules
# and do the work
include ${K1_TOOLCHAIN_DIR}/share/make/Makefile.mppaipc


# Define your rules after the include, otherwise, the first rule will be treated
# as the default rule
my_help_rule:
	@echo Additional Targets:
	@echo "  clean-traces    - clean generated traces"
	@echo "  run-sim         - run application on simulator (using only one MPPA)"
	@echo "  run-hw          - run application on board"
	@echo "                    when local, it uses config from /mppa/board*, assumes all boards are the same" 
	@echo "                    when remote, it needs REMOTE_ADDRESS/REMOTE_USER to be set" 

post_build_rule:
	@echo "Build done"

.PHONY: run-sim
run-sim: all
	${K1_TOOLCHAIN_DIR}/bin/k1-pciesim-runner -a ${arch} ${BIN_DIR}/${HOST_APPLI_NAME}  -b 1 -m 1 -i 1 ${BIN_DIR}/${multibin-name} ${K1_IO_APPLI_NAME} ${K1_CLUSTER_APPLI_NAME}

.PHONY: run-hw
ifeq (${remote}, true)
# When running on a remote target, the following variable must be defined.
# k1-remote-runner will also need a password. We STRONGLY recommend
# to do not use --passwd option for that.
# Instead, please consider providing the proper user ssh public key
# from developer to the remote target (you may use ssh-copyid standard tool).
# It will allow connecting without providing a password.
REMOTE_ADDRESS := $(shell echo $$REMOTE_ADDRESS)
REMOTE_USER := $(shell echo $$REMOTE_USER)
run-hw: all
	@[ "" != "${REMOTE_ADDRESS}" ]  || ( echo "REMOTE_ADDRESS not defined"  ; false)
	@[ "" != "${REMOTE_USER}" ]     || ( echo "REMOTE_USER not defined"     ; false)
	k1-remote-runner --remote ${REMOTE_ADDRESS} --user ${REMOTE_USER} --host_exec ${BIN_DIR}/${host-bin} --multibinary ${BIN_DIR}/${multibin-name} -- ${BIN_DIR}/${host-bin} -b 1 -m 1 -i 1 ${BIN_DIR}/${multibin-name} ${K1_IO_APPLI_NAME} ${K1_CLUSTER_APPLI_NAME}
else
# run-hw call the script "project_mppa.sh", so it auto-adapts to the available boards in /mppa/board*
run-hw: all
	./project_mppa.sh

run-hw-trace-no-visu: all
	${K1_TOOLCHAIN_DIR}/bin/k1-trace-util -a -l 0XFF -- ./project_mppa.sh
	if [ ! -s trace.dump.0 ] ; then echo "Error: no trace generated" ; false ; fi

endif

.PHONY: clean-traces
clean-traces:
	rm -rf MPPA.*.DSU.* andey.*.trace* bostan.*.trace* trace.dump.*
