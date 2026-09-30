# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

# parallel processing default
#VERSION = PARALLEL
# uncomment the next line to make serial processing the default
#VERSION = SERIAL

CXX = g++

CHECK_OPENMP_SRC = checkOpenMP/checkOpenMP.cpp
CHECK_OPENMP_OBJ = $(CHECK_OPENMP_SRC:cpp=o)
VERSION := $(shell $(CXX) -fopenmp -c $(CHECK_OPENMP_SRC) 2>/dev/null && ([ $$? -eq 0 ] && echo "PARALLEL") || echo "SERIAL")

ifeq "$(VERSION)" "PARALLEL"
CXXFLAGS = -Wall -fopenmp -O3
LINKER = -lm  -lstdc++ -fopenmp
else ifeq "$(VERSION)" "SERIAL"
CXXFLAGS = -Wall -O3
LINKER = -lm  -lstdc++
else #"$(VERSION)" "DEBUG"
CXXFLAGS = -Wall -O0 -g
LINKER = -g -lm  -lstdc++
endif

ENSCRIPT_SWITCHES = -2rjC -Ecpp --color --mark-wrapped-lines=arrow --line-numbers -H2
A2PS_SWITCHES = -r --columns=2 -j -C --tabsize=4 --pretty-print=c++ -g -b

#-------------------------------------
# Targets
#-------------------------------------
EXE = merced

OBJS = Compton.o Ecm_Elab_geom.o Energy_groups.o \
  Legendre2Body.o Legendre_base.o Legendre_data.o MadlandNix.o \
  Maxwell.o Vcm_Vlab_Hit.o Watt.o adapt_quad.o adapt_quad_info.o \
  angle_dist.o box_geom.o cm_Legendre.o \
  cm_doubleDiff.o coef_vector.o \
  coherent.o cumulative_points.o data_parser.o dd_vector.o \
  doubleDiff_base.o energy_dist.o \
  energy_dist_base.o energy_function.o evaporation.o expn.o \
  general_evap.o global_params.o igam.o isotropic.o \
  joint_dist.o joint_dist_data.o kalbach.o kalbach_data.o \
  lab_doubleDiff.o mappings.o \
  math_util.o messaging.o mtherr.o param_base.o phase_space.o \
  quad_methods.o reaction.o relativistic.o standard_Legendre.o \
  transfer.o two_step.o \
  two_step_hit.o uncorrelated.o x_vector.o main.o

SOURCES = main.cpp \
  Compton.cpp			Compton.hpp \
  Ecm_Elab_geom.cpp		Ecm_Elab_geom.hpp \
  Energy_groups.cpp 		Energy_groups.hpp \
  Legendre2Body.cpp		Legendre2Body.hpp \
  Legendre_base.cpp		Legendre_base.hpp \
  Legendre_data.cpp		Legendre_data.hpp \
  MadlandNix.cpp 		MadlandNix.hpp \
  Maxwell.cpp			Maxwell.hpp \
  Vcm_Vlab_Hit.cpp		Vcm_Vlab_Hit.hpp \
  Watt.cpp			Watt.hpp \
  adapt_quad.cpp		adapt_quad.hpp \
  adapt_quad_info.cpp 		adapt_quad_info.hpp \
  angle_dist.cpp		angle_dist.hpp \
  box_geom.cpp 			box_geom.hpp \
  cm_Legendre.cpp		cm_Legendre.hpp \
  cm_doubleDiff.cpp		cm_doubleDiff.hpp \
  coef_vector.cpp 		coef_vector.hpp \
  coherent.cpp			coherent.hpp \
  cumulative_points.cpp		cumulative_points.hpp \
  data_parser.cpp		data_parser.hpp \
  dd_vector.cpp			dd_vector.hpp \
  doubleDiff_base.cpp		doubleDiff_base.hpp \
  energy_dist.cpp 		energy_dist.hpp \
  energy_dist_base.cpp		energy_dist_base.hpp \
  energy_function.cpp		energy_function.hpp \
  evaporation.cpp		evaporation.hpp \
  expn.cpp \
  general_evap.cpp		general_evap.hpp \
  global_params.cpp		global_params.hpp \
  igam.cpp \
  isotropic.cpp 		isotropic.hpp \
  joint_dist.cpp		joint_dist.hpp \
  joint_dist_data.cpp		joint_dist_data.hpp \
  kalbach.cpp			kalbach.hpp \
  kalbach_data.cpp		kalbach_data.hpp \
  lab_doubleDiff.cpp		lab_doubleDiff.hpp \
  mappings.cpp 			mappings.hpp \
  math_util.cpp			math_util.hpp \
    max_output_order.hpp \
  messaging.cpp			messaging.hpp \
  mtherr.cpp \
  param_base.cpp		param_base.hpp \
  phase_space.cpp 		phase_space.hpp \
    protos.h  \
  quad_methods.cpp  		quad_methods.hpp \
  reaction.cpp			reaction.hpp \
  relativistic.cpp		relativistic.hpp \
  standard_Legendre.cpp		standard_Legendre.hpp \
  transfer.cpp			transfer.hpp \
  two_step.cpp			two_step.hpp \
  two_step_hit.cpp		two_step_hit.hpp \
  uncorrelated.cpp 		uncorrelated.hpp \
  x_vector.cpp			x_vector.hpp

#-------------------------------------
# Action item(s)
#-------------------------------------
all: exe 

default : all

exe: $(EXE)

clean:
	@$(RM) -f *.o $(CHECK_OPENMP_OBJ)

cleaner: clean
	@$(RM) -f *.*~ mf*mt* yo*c*i*s* 

realclean: cleaner
	@$(RM) -f $(EXE) $(LIBRARY)

install: all
	cp $(EXE) ../bin

uninstall:

doc:

help:

#-------------------------------------
#  Suffix rules 
#-------------------------------------
.SUFFIXES : .hpp_ps .cpp_ps .cpp .hpp .c .h .o 

.hpp.hpp_ps :
#	 a2ps $(A2PS_SWITCHES) -o $@ $<
	 sed 's/	/  /g' $< | enscript $(ENSCRIPT_SWITCHES) -o $@
.cpp.cpp_ps :
#	 a2ps $(A2PS_SWITCHES) -o $@ $<
	 sed 's/	/  /g' $< | enscript $(ENSCRIPT_SWITCHES) -o $@
.c.o :
	$(CXX) $(CXXFLAGS) -c -o $@ $< 
.cpp.o :
	$(CXX) $(CXXFLAGS) -c -o $@ $< 
.o :
	$(CXX) -o $@ $<

#-------------------------------------
# Rules
#-------------------------------------

$(EXE): $(OBJS)
	$(CXX) $(OBJS) $(LINKER) -o $(EXE)

version.hpp: $(SOURCES)
	python create_version_hpp.py

Compton.o: Compton.cpp \
  Compton.hpp Energy_groups.hpp Legendre_base.hpp \
  Legendre_data.hpp adapt_quad.hpp adapt_quad_info.hpp \
  coef_vector.hpp data_parser.hpp dd_vector.hpp \
  global_params.hpp math_util.hpp max_output_order.hpp \
  messaging.hpp param_base.hpp quad_methods.hpp quad_param.hpp \
  transfer.hpp x_vector.hpp

Ecm_Elab_geom.o: Ecm_Elab_geom.cpp \
  Ecm_Elab_geom.hpp Energy_groups.hpp Legendre_base.hpp \
  Legendre_data.hpp Vcm_Vlab_Hit.hpp box_geom.hpp \
  coef_vector.hpp data_parser.hpp dd_vector.hpp \
  global_params.hpp mappings.hpp math_util.hpp \
  max_output_order.hpp messaging.hpp param_base.hpp \
  quad_param.hpp

Energy_groups.o: Energy_groups.cpp \
  Energy_groups.hpp data_parser.hpp global_params.hpp \
  messaging.hpp

Legendre2Body.o: Legendre2Body.cpp \
  Energy_groups.hpp Legendre2Body.hpp Legendre_base.hpp \
  Legendre_data.hpp adapt_quad.hpp adapt_quad_info.hpp \
  angle_dist.hpp box_geom.hpp coef_vector.hpp data_parser.hpp \
  dd_vector.hpp global_params.hpp mappings.hpp math_util.hpp \
  max_output_order.hpp messaging.hpp param_base.hpp \
  quad_methods.hpp quad_param.hpp relativistic.hpp transfer.hpp

Legendre_base.o: Legendre_base.cpp \
  Legendre_base.hpp messaging.hpp

Legendre_data.o: Legendre_data.cpp \
  Energy_groups.hpp Legendre_base.hpp Legendre_data.hpp \
  adapt_quad.hpp adapt_quad_info.hpp coef_vector.hpp \
  data_parser.hpp dd_vector.hpp global_params.hpp math_util.hpp \
  max_output_order.hpp messaging.hpp param_base.hpp \
  quad_methods.hpp quad_param.hpp

MadlandNix.o: MadlandNix.cpp \
  Energy_groups.hpp Legendre_base.hpp Legendre_data.hpp \
  MadlandNix.hpp box_geom.hpp coef_vector.hpp data_parser.hpp \
  dd_vector.hpp energy_function.hpp global_params.hpp \
  math_util.hpp max_output_order.hpp messaging.hpp \
  param_base.hpp protos.h quad_methods.hpp quad_param.hpp \
  transfer.hpp

Maxwell.o: Maxwell.cpp \
  Energy_groups.hpp Legendre_base.hpp Legendre_data.hpp \
  Maxwell.hpp adapt_quad.hpp adapt_quad_info.hpp box_geom.hpp \
  coef_vector.hpp data_parser.hpp dd_vector.hpp \
  energy_function.hpp evaporation.hpp global_params.hpp \
  max_output_order.hpp messaging.hpp param_base.hpp protos.h \
  quad_methods.hpp quad_param.hpp transfer.hpp

Vcm_Vlab_Hit.o: Vcm_Vlab_Hit.cpp \
  Energy_groups.hpp Legendre_base.hpp Legendre_data.hpp \
  Vcm_Vlab_Hit.hpp box_geom.hpp coef_vector.hpp data_parser.hpp \
  dd_vector.hpp global_params.hpp mappings.hpp math_util.hpp \
  max_output_order.hpp messaging.hpp param_base.hpp \
  quad_param.hpp

Watt.o: Watt.cpp \
  Energy_groups.hpp Legendre_base.hpp Legendre_data.hpp \
  Watt.hpp adapt_quad.hpp adapt_quad_info.hpp box_geom.hpp \
  coef_vector.hpp data_parser.hpp dd_vector.hpp \
  energy_function.hpp global_params.hpp max_output_order.hpp \
  messaging.hpp param_base.hpp protos.h quad_methods.hpp \
  quad_param.hpp transfer.hpp

adapt_quad.o: adapt_quad.cpp \
  Energy_groups.hpp Legendre_base.hpp Legendre_data.hpp \
  adapt_quad.hpp adapt_quad_info.hpp coef_vector.hpp \
  data_parser.hpp dd_vector.hpp global_params.hpp \
  max_output_order.hpp messaging.hpp param_base.hpp \
  quad_methods.hpp quad_param.hpp

adapt_quad_info.o: adapt_quad_info.cpp \
  Legendre_base.hpp adapt_quad_info.hpp coef_vector.hpp \
  global_params.hpp max_output_order.hpp

angle_dist.o: angle_dist.cpp \
  Energy_groups.hpp Legendre_base.hpp Legendre_data.hpp \
  adapt_quad.hpp adapt_quad_info.hpp angle_dist.hpp \
  box_geom.hpp coef_vector.hpp data_parser.hpp dd_vector.hpp \
  global_params.hpp mappings.hpp math_util.hpp \
  max_output_order.hpp messaging.hpp param_base.hpp \
  quad_methods.hpp quad_param.hpp relativistic.hpp transfer.hpp

box_geom.o: box_geom.cpp \
  box_geom.hpp data_parser.hpp dd_vector.hpp global_params.hpp \
  messaging.hpp

cm_Legendre.o: cm_Legendre.cpp \
  Ecm_Elab_geom.hpp Energy_groups.hpp Legendre_base.hpp \
  Legendre_data.hpp Vcm_Vlab_Hit.hpp adapt_quad.hpp \
  adapt_quad_info.hpp box_geom.hpp cm_Legendre.hpp \
  coef_vector.hpp cumulative_points.hpp data_parser.hpp \
  dd_vector.hpp global_params.hpp mappings.hpp math_util.hpp \
  max_output_order.hpp messaging.hpp param_base.hpp \
  quad_methods.hpp quad_param.hpp transfer.hpp

cm_doubleDiff.o: cm_doubleDiff.cpp \
  Ecm_Elab_geom.hpp Energy_groups.hpp Legendre_base.hpp \
  Legendre_data.hpp Vcm_Vlab_Hit.hpp adapt_quad.hpp \
  adapt_quad_info.hpp box_geom.hpp cm_doubleDiff.hpp \
  coef_vector.hpp data_parser.hpp dd_vector.hpp \
  doubleDiff_base.hpp global_params.hpp mappings.hpp \
  math_util.hpp max_output_order.hpp messaging.hpp \
  param_base.hpp quad_methods.hpp quad_param.hpp transfer.hpp

coef_vector.o: coef_vector.cpp \
  Legendre_base.hpp coef_vector.hpp global_params.hpp \
  max_output_order.hpp messaging.hpp

coherent.o: coherent.cpp \
  Energy_groups.hpp Legendre_base.hpp Legendre_data.hpp \
  adapt_quad.hpp adapt_quad_info.hpp box_geom.hpp \
  coef_vector.hpp coherent.hpp data_parser.hpp dd_vector.hpp \
  global_params.hpp math_util.hpp max_output_order.hpp \
  messaging.hpp param_base.hpp quad_methods.hpp quad_param.hpp \
  transfer.hpp x_vector.hpp

cumulative_points.o: cumulative_points.cpp \
  cumulative_points.hpp data_parser.hpp dd_vector.hpp \
  global_params.hpp messaging.hpp

data_parser.o: data_parser.cpp \
  data_parser.hpp messaging.hpp

dd_vector.o: dd_vector.cpp \
  data_parser.hpp dd_vector.hpp global_params.hpp messaging.hpp

doubleDiff_base.o: doubleDiff_base.cpp \
  data_parser.hpp dd_vector.hpp doubleDiff_base.hpp \
  global_params.hpp messaging.hpp

energy_dist.o: energy_dist.cpp \
  Energy_groups.hpp Legendre_base.hpp Legendre_data.hpp \
  box_geom.hpp coef_vector.hpp cumulative_points.hpp \
  data_parser.hpp dd_vector.hpp energy_dist.hpp \
  energy_dist_base.hpp global_params.hpp math_util.hpp \
  max_output_order.hpp messaging.hpp param_base.hpp \
  quad_methods.hpp quad_param.hpp standard_Legendre.hpp \
  transfer.hpp

energy_dist_base.o: energy_dist_base.cpp \
  Energy_groups.hpp Legendre_base.hpp Legendre_data.hpp \
  box_geom.hpp coef_vector.hpp cumulative_points.hpp \
  data_parser.hpp dd_vector.hpp energy_dist_base.hpp \
  global_params.hpp max_output_order.hpp messaging.hpp \
  param_base.hpp quad_methods.hpp quad_param.hpp transfer.hpp

energy_function.o: energy_function.cpp \
  Energy_groups.hpp Legendre_base.hpp Legendre_data.hpp \
  adapt_quad.hpp adapt_quad_info.hpp box_geom.hpp \
  coef_vector.hpp data_parser.hpp dd_vector.hpp \
  energy_function.hpp global_params.hpp max_output_order.hpp \
  messaging.hpp param_base.hpp quad_methods.hpp quad_param.hpp \
  transfer.hpp

evaporation.o: evaporation.cpp \
  Energy_groups.hpp Legendre_base.hpp Legendre_data.hpp \
  adapt_quad.hpp adapt_quad_info.hpp box_geom.hpp \
  coef_vector.hpp data_parser.hpp dd_vector.hpp \
  energy_function.hpp evaporation.hpp global_params.hpp \
  max_output_order.hpp messaging.hpp param_base.hpp protos.h \
  quad_methods.hpp quad_param.hpp transfer.hpp

expn.o: expn.cpp protos.h

general_evap.o: general_evap.cpp \
  Energy_groups.hpp Legendre_base.hpp Legendre_data.hpp \
  adapt_quad.hpp adapt_quad_info.hpp coef_vector.hpp \
  data_parser.hpp dd_vector.hpp general_evap.hpp \
  global_params.hpp max_output_order.hpp messaging.hpp \
  param_base.hpp quad_methods.hpp quad_param.hpp transfer.hpp

global_params.o: global_params.cpp \
  data_parser.hpp global_params.hpp messaging.hpp

igam.o: igam.cpp protos.h

isotropic.o: isotropic.cpp \
  Energy_groups.hpp Legendre_base.hpp Legendre_data.hpp \
  box_geom.hpp coef_vector.hpp cumulative_points.hpp \
  data_parser.hpp dd_vector.hpp energy_dist.hpp \
  energy_dist_base.hpp global_params.hpp isotropic.hpp \
  math_util.hpp max_output_order.hpp messaging.hpp \
  param_base.hpp quad_methods.hpp quad_param.hpp \
  standard_Legendre.hpp transfer.hpp

joint_dist.o: joint_dist.cpp \
  Energy_groups.hpp Legendre_base.hpp Legendre_data.hpp \
  adapt_quad.hpp adapt_quad_info.hpp angle_dist.hpp \
  box_geom.hpp coef_vector.hpp cumulative_points.hpp \
  data_parser.hpp dd_vector.hpp global_params.hpp \
  joint_dist.hpp joint_dist_data.hpp mappings.hpp math_util.hpp \
  max_output_order.hpp messaging.hpp param_base.hpp \
  quad_methods.hpp quad_param.hpp relativistic.hpp transfer.hpp

joint_dist_data.o: joint_dist_data.cpp \
  data_parser.hpp dd_vector.hpp global_params.hpp \
  joint_dist_data.hpp messaging.hpp

kalbach.o: kalbach.cpp \
  Ecm_Elab_geom.hpp Energy_groups.hpp Legendre_base.hpp \
  Legendre_data.hpp Vcm_Vlab_Hit.hpp adapt_quad.hpp \
  adapt_quad_info.hpp box_geom.hpp coef_vector.hpp \
  data_parser.hpp dd_vector.hpp global_params.hpp kalbach.hpp \
  kalbach_data.hpp mappings.hpp math_util.hpp \
  max_output_order.hpp messaging.hpp param_base.hpp \
  quad_methods.hpp quad_param.hpp transfer.hpp

kalbach_data.o: kalbach_data.cpp \
  data_parser.hpp dd_vector.hpp global_params.hpp \
  kalbach_data.hpp mappings.hpp messaging.hpp

lab_doubleDiff.o: lab_doubleDiff.cpp \
  Energy_groups.hpp Legendre_base.hpp Legendre_data.hpp \
  box_geom.hpp coef_vector.hpp cumulative_points.hpp \
  data_parser.hpp dd_vector.hpp doubleDiff_base.hpp \
  global_params.hpp lab_doubleDiff.hpp math_util.hpp \
  max_output_order.hpp param_base.hpp quad_methods.hpp \
  quad_param.hpp standard_Legendre.hpp transfer.hpp

mappings.o: mappings.cpp \
  Legendre_data.hpp coef_vector.hpp data_parser.hpp \
  dd_vector.hpp global_params.hpp mappings.hpp math_util.hpp \
  max_output_order.hpp messaging.hpp

math_util.o: math_util.cpp \
  Legendre_data.hpp coef_vector.hpp data_parser.hpp \
  dd_vector.hpp global_params.hpp math_util.hpp \
  max_output_order.hpp messaging.hpp

messaging.o: messaging.cpp \
  global_params.hpp messaging.hpp

mtherr.o: mtherr.cpp protos.h

param_base.o: param_base.cpp \
  Energy_groups.hpp Legendre_base.hpp Legendre_data.hpp \
  coef_vector.hpp data_parser.hpp dd_vector.hpp \
  global_params.hpp max_output_order.hpp messaging.hpp \
  param_base.hpp quad_param.hpp

phase_space.o: phase_space.cpp \
  Ecm_Elab_geom.hpp Energy_groups.hpp Legendre_base.hpp \
  Legendre_data.hpp Vcm_Vlab_Hit.hpp adapt_quad.hpp \
  adapt_quad_info.hpp box_geom.hpp coef_vector.hpp \
  data_parser.hpp dd_vector.hpp global_params.hpp mappings.hpp \
  math_util.hpp max_output_order.hpp messaging.hpp \
  param_base.hpp phase_space.hpp quad_methods.hpp \
  quad_param.hpp transfer.hpp

quad_methods.o: quad_methods.cpp \
  Energy_groups.hpp Legendre_base.hpp Legendre_data.hpp \
  coef_vector.hpp data_parser.hpp dd_vector.hpp \
  max_output_order.hpp param_base.hpp quad_methods.hpp \
  quad_param.hpp

reaction.o: reaction.cpp \
  Compton.hpp Ecm_Elab_geom.hpp Energy_groups.hpp \
  Legendre2Body.hpp Legendre_base.hpp Legendre_data.hpp \
  MadlandNix.hpp Maxwell.hpp Vcm_Vlab_Hit.hpp Watt.hpp \
  angle_dist.hpp box_geom.hpp cm_Legendre.hpp cm_doubleDiff.hpp \
  coef_vector.hpp coherent.hpp cumulative_points.hpp \
  data_parser.hpp dd_vector.hpp doubleDiff_base.hpp \
  energy_dist.hpp energy_dist_base.hpp energy_function.hpp \
  evaporation.hpp general_evap.hpp global_params.hpp \
  isotropic.hpp joint_dist.hpp joint_dist_data.hpp kalbach.hpp \
  kalbach_data.hpp lab_doubleDiff.hpp mappings.hpp \
  math_util.hpp max_output_order.hpp messaging.hpp \
  param_base.hpp phase_space.hpp quad_methods.hpp \
  quad_param.hpp reaction.hpp relativistic.hpp \
  standard_Legendre.hpp transfer.hpp \
  two_step.hpp two_step_hit.hpp uncorrelated.hpp x_vector.hpp

relativistic.o: relativistic.cpp \
  Legendre_data.hpp coef_vector.hpp data_parser.hpp \
  dd_vector.hpp global_params.hpp mappings.hpp math_util.hpp \
  max_output_order.hpp messaging.hpp relativistic.hpp

standard_Legendre.o: standard_Legendre.cpp \
  Energy_groups.hpp Legendre_base.hpp Legendre_data.hpp \
  adapt_quad.hpp adapt_quad_info.hpp box_geom.hpp \
  coef_vector.hpp cumulative_points.hpp data_parser.hpp \
  dd_vector.hpp global_params.hpp math_util.hpp \
  max_output_order.hpp messaging.hpp param_base.hpp \
  quad_methods.hpp quad_param.hpp standard_Legendre.hpp \
  transfer.hpp

transfer.o: transfer.cpp \
  Energy_groups.hpp Legendre_base.hpp Legendre_data.hpp \
  coef_vector.hpp data_parser.hpp dd_vector.hpp \
  global_params.hpp max_output_order.hpp messaging.hpp \
  param_base.hpp quad_methods.hpp quad_param.hpp transfer.hpp

two_step.o: two_step.cpp \
  Energy_groups.hpp Legendre_base.hpp Legendre_data.hpp \
  adapt_quad.hpp adapt_quad_info.hpp box_geom.hpp \
  coef_vector.hpp data_parser.hpp dd_vector.hpp \
  global_params.hpp mappings.hpp math_util.hpp \
  max_output_order.hpp messaging.hpp param_base.hpp \
  quad_methods.hpp quad_param.hpp relativistic.hpp transfer.hpp \
  two_step.hpp two_step_hit.hpp

two_step_hit.o: two_step_hit.cpp \
  box_geom.hpp data_parser.hpp dd_vector.hpp global_params.hpp \
  mappings.hpp messaging.hpp relativistic.hpp two_step_hit.hpp

uncorrelated.o: uncorrelated.cpp \
  Energy_groups.hpp Legendre2Body.hpp Legendre_base.hpp \
  Legendre_data.hpp adapt_quad.hpp adapt_quad_info.hpp \
  angle_dist.hpp box_geom.hpp coef_vector.hpp \
  cumulative_points.hpp data_parser.hpp dd_vector.hpp \
  energy_dist.hpp energy_dist_base.hpp global_params.hpp \
  mappings.hpp math_util.hpp max_output_order.hpp messaging.hpp \
  param_base.hpp quad_methods.hpp quad_param.hpp \
  relativistic.hpp standard_Legendre.hpp transfer.hpp \
  uncorrelated.hpp

x_vector.o: x_vector.cpp \
  data_parser.hpp dd_vector.hpp global_params.hpp messaging.hpp \
  x_vector.hpp

main.o: main.cpp \
  Compton.hpp Ecm_Elab_geom.hpp Energy_groups.hpp \
  Legendre2Body.hpp Legendre_base.hpp Legendre_data.hpp \
  MadlandNix.hpp Maxwell.hpp Vcm_Vlab_Hit.hpp Watt.hpp \
  angle_dist.hpp box_geom.hpp cm_Legendre.hpp cm_doubleDiff.hpp \
  coef_vector.hpp coherent.hpp cumulative_points.hpp \
  data_parser.hpp dd_vector.hpp doubleDiff_base.hpp \
  energy_dist.hpp energy_dist_base.hpp energy_function.hpp \
  evaporation.hpp general_evap.hpp global_params.hpp \
  isotropic.hpp joint_dist.hpp joint_dist_data.hpp kalbach.hpp \
  kalbach_data.hpp lab_doubleDiff.hpp mappings.hpp \
  math_util.hpp max_output_order.hpp messaging.hpp \
  param_base.hpp phase_space.hpp quad_methods.hpp \
  quad_param.hpp reaction.hpp relativistic.hpp \
  standard_Legendre.hpp transfer.hpp \
  two_step.hpp two_step_hit.hpp uncorrelated.hpp x_vector.hpp
