include ../machine/pluton.makefile
.SUFFIXES:
.SUFFIXES: .f90 .o
.f90.o:
	$(FC) -c $(FFLAGS) $<
.SUFFIXES: .F90 .o
.F90.o:
	$(FC) -c $(FFLAGS) $<

wrapper: compil_verbose

compil_verbose: elasto.exe
	@echo -e "\033[1;31m===== Compilation successfull =====\033[0m"
	@cd ../.. ;ln -fs builddir/elasto.exe ../elasto.exe ; rm ELASTO/builddir/*.mod; rm ELASTO/builddir/*.o
	@echo -e "\033[1;30mProgram is in `(cd ../..;pwd)`/scaleExe\033[0m"
FOBJ=avg_expectation_tools.o avs_datalayout.o cart_mesh_tools.o cart_topology.o data.o datalayout.o dealias_tools.o energy_transfers.o energy_transfers_lib.o EqDM_lorentz_models.o fileio.o filtering_tools.o forcing.o forcing_scalar.o geometric_shape_tools.o interface_layout_mesh.o plasticflow.o maindatalayout.o main.o mpilayout_tools.o parallel_tools.o param.o parser_tools.o physical_values_tools.o post_HD_ener_lib.o post_hd.o post_HD_kinHel_lib.o post_hd_out.o post_in_simu.o post_kolmo.o post_lib2D.o post_lib.o post_MHD_ener_lib.o post_mhd.o post_MHD_helic_lib.o post_mhd_out.o post_specific_to_setup.o post_visus.o precision_tools.o random_tools.o rediscretization_tools.o solver.o sort_find_tools.o stat_tools.o structure_tools.o surcharge_vtkxml.o system_tools.o transforms_tools.o vel_models.o vtkxml_bin.o vtkxml.o wavenumber_tools.o avs_complex.o avs_real.o communicators_tools.o complexdatalayout.o complexparallel_tools.o complexVector_datalayout.o complexVector_parallel_tools.o differential_tools.o hdf5.o hdf5_wrapper_cartmesh.o hdf5_wrapper_datalayout.o interface_fann.o interface_io.o interpolation_velo.o mhd_model.o options.o post_out_simu.o post_scalar.o post_velocity.o realdatalayout.o realparallel_tools.o realVector_datalayout.o realVector_parallel_tools.o sca_models.o sgs_model_tools.o toolbox.o xdmf.o fann_run.o

elasto.exe: $(FOBJ)
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $(FOBJ) $(LIBS)

avg_expectation_tools.o : avg_expectation_tools.f90 toolbox.o sort_find_tools.o mpilayout_tools.o stat_tools.o precision_tools.o parallel_tools.o datalayout.o 
avs_datalayout.o : avs_datalayout.f90 fileio.o parallel_tools.o datalayout.o avs_complex.o avs_real.o precision_tools.o 
cart_mesh_tools.o : cart_mesh_tools.f90 precision_tools.o 
cart_topology.o : cart_topology.f90 cart_mesh_tools.o precision_tools.o 
data.o : data.f90 interface_io.o parallel_tools.o parser_tools.o cart_topology.o param.o wavenumber_tools.o datalayout.o 
datalayout.o : datalayout.f90 realVector_datalayout.o complexVector_datalayout.o realdatalayout.o complexdatalayout.o 
dealias_tools.o : dealias_tools.f90 data.o wavenumber_tools.o datalayout.o precision_tools.o 
energy_transfers.o : energy_transfers.f90 energy_transfers_lib.o transforms_tools.o data.o fileio.o parallel_tools.o datalayout.o wavenumber_tools.o 
energy_transfers_lib.o : energy_transfers_lib.f90 fileio.o parallel_tools.o datalayout.o wavenumber_tools.o 
EqDM_lorentz_models.o : EqDM_lorentz_models.f90 differential_tools.o datalayout.o sgs_model_tools.o post_lib.o toolbox.o transforms_tools.o wavenumber_tools.o 
fileio.o : fileio.f90 
filtering_tools.o : filtering_tools.f90 datalayout.o parallel_tools.o stat_tools.o transforms_tools.o wavenumber_tools.o precision_tools.o 
forcing.o : forcing.f90 random_tools.o transforms_tools.o data.o communicators_tools.o parallel_tools.o param.o datalayout.o wavenumber_tools.o precision_tools.o 
forcing_scalar.o : forcing_scalar.f90 random_tools.o datalayout.o wavenumber_tools.o precision_tools.o 
geometric_shape_tools.o : geometric_shape_tools.f90 parallel_tools.o physical_values_tools.o stat_tools.o realdatalayout.o precision_tools.o 
interface_layout_mesh.o : interface_layout_mesh.f90 realdatalayout.o cart_mesh_tools.o precision_tools.o 
plasticflow.o : plasticflow.f90 precision_tools.o datalayout.o wavenumber_tools.o
maindatalayout.o : maindatalayout.f90 mpilayout_tools.o precision_tools.o 
main.o : main.f90 post_out_simu.o post_in_simu.o post_lib.o solver.o avs_datalayout.o precision_tools.o datalayout.o data.o random_tools.o communicators_tools.o mpilayout_tools.o param.o options.o 
mpilayout_tools.o : mpilayout_tools.f90 options.o 
parallel_tools.o : parallel_tools.f90 mpilayout_tools.o communicators_tools.o precision_tools.o realVector_parallel_tools.o complexVector_parallel_tools.o realparallel_tools.o complexparallel_tools.o 
param.o : param.f90 parser_tools.o precision_tools.o 
parser_tools.o : parser_tools.f90 fileio.o precision_tools.o 
physical_values_tools.o : physical_values_tools.f90 avs_datalayout.o filtering_tools.o differential_tools.o wavenumber_tools.o datalayout.o stat_tools.o precision_tools.o 
post_HD_ener_lib.o : post_HD_ener_lib.f90 physical_values_tools.o vel_models.o mhd_model.o filtering_tools.o sgs_model_tools.o transforms_tools.o wavenumber_tools.o differential_tools.o stat_tools.o toolbox.o datalayout.o 
post_hd.o : post_hd.f90 post_HD_ener_lib.o vel_models.o mhd_model.o filtering_tools.o sgs_model_tools.o transforms_tools.o wavenumber_tools.o differential_tools.o stat_tools.o toolbox.o datalayout.o 
post_HD_kinHel_lib.o : post_HD_kinHel_lib.f90 physical_values_tools.o vel_models.o mhd_model.o filtering_tools.o sgs_model_tools.o transforms_tools.o wavenumber_tools.o differential_tools.o stat_tools.o toolbox.o datalayout.o 
post_hd_out.o : post_hd_out.f90 data.o post_HD_kinHel_lib.o post_HD_ener_lib.o vel_models.o filtering_tools.o sgs_model_tools.o transforms_tools.o wavenumber_tools.o differential_tools.o stat_tools.o toolbox.o datalayout.o 
post_in_simu.o : post_in_simu.f90 interface_io.o filtering_tools.o toolbox.o rediscretization_tools.o differential_tools.o post_velocity.o cart_topology.o data.o param.o post_hd.o parser_tools.o stat_tools.o system_tools.o post_kolmo.o post_mhd.o surcharge_vtkxml.o avs_datalayout.o post_lib2D.o post_scalar.o post_lib.o datalayout.o precision_tools.o 
post_kolmo.o : post_kolmo.f90 fileio.o post_lib2D.o transforms_tools.o wavenumber_tools.o differential_tools.o stat_tools.o toolbox.o datalayout.o 
post_lib2D.o : post_lib2D.f90 fileio.o parallel_tools.o data.o transforms_tools.o stat_tools.o toolbox.o avs_datalayout.o datalayout.o precision_tools.o 
post_lib.o : post_lib.f90 mpilayout_tools.o interface_io.o physical_values_tools.o precision_tools.o parser_tools.o stat_tools.o toolbox.o rediscretization_tools.o data.o surcharge_vtkxml.o sgs_model_tools.o avs_datalayout.o wavenumber_tools.o transforms_tools.o datalayout.o parallel_tools.o fileio.o 
post_MHD_ener_lib.o : post_MHD_ener_lib.f90 vel_models.o mhd_model.o filtering_tools.o sgs_model_tools.o transforms_tools.o wavenumber_tools.o differential_tools.o stat_tools.o toolbox.o datalayout.o 
post_mhd.o : post_mhd.f90 EqDM_lorentz_models.o vel_models.o mhd_model.o filtering_tools.o sgs_model_tools.o transforms_tools.o wavenumber_tools.o differential_tools.o stat_tools.o toolbox.o datalayout.o 
post_MHD_helic_lib.o : post_MHD_helic_lib.f90 vel_models.o mhd_model.o filtering_tools.o sgs_model_tools.o transforms_tools.o wavenumber_tools.o differential_tools.o stat_tools.o toolbox.o datalayout.o 
post_mhd_out.o : post_mhd_out.f90 post_MHD_ener_lib.o vel_models.o mhd_model.o filtering_tools.o sgs_model_tools.o transforms_tools.o wavenumber_tools.o differential_tools.o stat_tools.o toolbox.o datalayout.o data.o 
post_specific_to_setup.o : post_specific_to_setup.f90 transforms_tools.o post_scalar.o post_velocity.o post_lib.o datalayout.o differential_tools.o stat_tools.o toolbox.o parser_tools.o precision_tools.o 
post_visus.o : post_visus.f90 filtering_tools.o sgs_model_tools.o parallel_tools.o wavenumber_tools.o post_scalar.o interface_io.o mpilayout_tools.o post_lib.o transforms_tools.o physical_values_tools.o post_hd_out.o data.o differential_tools.o stat_tools.o toolbox.o datalayout.o avs_datalayout.o 
precision_tools.o : precision_tools.f90 
random_tools.o : random_tools.f90 precision_tools.o 
rediscretization_tools.o : rediscretization_tools.f90 parallel_tools.o mpilayout_tools.o transforms_tools.o datalayout.o 
solver.o : solver.f90 dealias_tools.o cart_topology.o interpolation_velo.o forcing_scalar.o forcing.o rediscretization_tools.o data.o parallel_tools.o wavenumber_tools.o plasticflow.o param.o datalayout.o precision_tools.o differential_tools.o 
sort_find_tools.o : sort_find_tools.f90 mpilayout_tools.o parallel_tools.o stat_tools.o datalayout.o precision_tools.o 
stat_tools.o : stat_tools.f90 toolbox.o fileio.o differential_tools.o transforms_tools.o wavenumber_tools.o mpilayout_tools.o datalayout.o parallel_tools.o precision_tools.o 
structure_tools.o : structure_tools.f90 precision_tools.o 
surcharge_vtkxml.o : surcharge_vtkxml.f90 interface_layout_mesh.o realdatalayout.o vtkxml.o precision_tools.o 
system_tools.o : system_tools.f90 
transforms_tools.o : transforms_tools.f90 datalayout.o parallel_tools.o 
vel_models.o : vel_models.f90 mpilayout_tools.o physical_values_tools.o filtering_tools.o sgs_model_tools.o wavenumber_tools.o transforms_tools.o differential_tools.o stat_tools.o toolbox.o datalayout.o 
vtkxml_bin.o : vtkxml_bin.f90 cart_topology.o cart_mesh_tools.o precision_tools.o 
vtkxml.o : vtkxml.f90 realdatalayout.o cart_topology.o cart_mesh_tools.o precision_tools.o 
wavenumber_tools.o : wavenumber_tools.f90 precision_tools.o 
avs_complex.o : avs_complex.F90 communicators_tools.o datalayout.o 
avs_real.o : avs_real.F90 communicators_tools.o datalayout.o 
communicators_tools.o : communicators_tools.F90 cart_topology.o precision_tools.o 
complexdatalayout.o : complexdatalayout.F90 parser_tools.o mpilayout_tools.o precision_tools.o maindatalayout.o 
complexparallel_tools.o : complexparallel_tools.F90 mpilayout_tools.o precision_tools.o communicators_tools.o datalayout.o 
complexVector_datalayout.o : complexVector_datalayout.F90 complexdatalayout.o parser_tools.o mpilayout_tools.o precision_tools.o maindatalayout.o 
complexVector_parallel_tools.o : complexVector_parallel_tools.F90 complexVector_datalayout.o communicators_tools.o datalayout.o complexparallel_tools.o 
differential_tools.o : differential_tools.F90 rediscretization_tools.o parallel_tools.o transforms_tools.o wavenumber_tools.o datalayout.o precision_tools.o 
hdf5.o : hdf5.F90 precision_tools.o 
hdf5_wrapper_cartmesh.o : hdf5_wrapper_cartmesh.F90 precision_tools.o cart_topology.o 
hdf5_wrapper_datalayout.o : hdf5_wrapper_datalayout.F90 datalayout.o hdf5.o precision_tools.o 
interface_fann.o : interface_fann.F90 datalayout.o 
interface_io.o : interface_io.F90 avs_datalayout.o rediscretization_tools.o mpilayout_tools.o datalayout.o hdf5_wrapper_datalayout.o avs_real.o communicators_tools.o hdf5.o parser_tools.o precision_tools.o 
interpolation_velo.o : interpolation_velo.F90 cart_topology.o precision_tools.o 
mhd_model.o : mhd_model.F90 differential_tools.o datalayout.o sgs_model_tools.o post_lib.o toolbox.o transforms_tools.o wavenumber_tools.o 
options.o : options.F90 
post_out_simu.o : post_out_simu.F90 system_tools.o interface_io.o energy_transfers.o post_mhd_out.o post_visus.o post_hd_out.o stat_tools.o post_specific_to_setup.o cart_topology.o surcharge_vtkxml.o wavenumber_tools.o post_velocity.o post_scalar.o post_lib.o data.o datalayout.o precision_tools.o param.o avg_expectation_tools.o 
post_scalar.o : post_scalar.F90 avs_datalayout.o mpilayout_tools.o geometric_shape_tools.o physical_values_tools.o filtering_tools.o differential_tools.o parser_tools.o fileio.o dealias_tools.o sgs_model_tools.o parallel_tools.o avg_expectation_tools.o wavenumber_tools.o sca_models.o post_lib.o transforms_tools.o stat_tools.o toolbox.o datalayout.o precision_tools.o 
post_velocity.o : post_velocity.F90 param.o data.o interpolation_velo.o sgs_model_tools.o avg_expectation_tools.o parallel_tools.o fileio.o differential_tools.o filtering_tools.o mpilayout_tools.o vel_models.o wavenumber_tools.o physical_values_tools.o post_lib.o transforms_tools.o stat_tools.o toolbox.o datalayout.o precision_tools.o 
realdatalayout.o : realdatalayout.F90 parser_tools.o mpilayout_tools.o precision_tools.o maindatalayout.o 
realparallel_tools.o : realparallel_tools.F90 mpilayout_tools.o precision_tools.o communicators_tools.o datalayout.o 
realVector_datalayout.o : realVector_datalayout.F90 realdatalayout.o parser_tools.o mpilayout_tools.o precision_tools.o maindatalayout.o 
realVector_parallel_tools.o : realVector_parallel_tools.F90 realVector_datalayout.o communicators_tools.o datalayout.o realparallel_tools.o 
sca_models.o : sca_models.F90 physical_values_tools.o filtering_tools.o avg_expectation_tools.o sgs_model_tools.o wavenumber_tools.o transforms_tools.o differential_tools.o stat_tools.o toolbox.o datalayout.o 
sgs_model_tools.o : sgs_model_tools.F90 mpilayout_tools.o physical_values_tools.o filtering_tools.o datalayout.o wavenumber_tools.o transforms_tools.o differential_tools.o stat_tools.o toolbox.o 
toolbox.o : toolbox.F90 wavenumber_tools.o mpilayout_tools.o datalayout.o parallel_tools.o precision_tools.o 
xdmf.o : xdmf.F90 fileio.o parser_tools.o precision_tools.o 
fann_run.o : fann_run.c
	$(CC) -c $(CFLAGS) fann_run.c

clean:
	@rm -rf *.o *.mod core ../../doc/doxygen/codescalar_doc *.exe ../*.exe
doc:
	@cd ../../doc/doxygen && doxygen Doxyfile && echo "Doc start in `pwd`codescalar_doc/html/index.html"
fortran_file:
	@sed -e "s/TAG_/REAL/g" -e "s/tag_/real/g"  < ../toolboxes/implementVector_parallel_tools.Fortran > `echo implementVector_parallel_tools.Fortran|sed -e "s/implement/real/" -e "s/Fortran/F90/"`
	@sed -e "s/TAG_/COMPLEX/g" -e "s/tag_/cmplx/g"  < ../toolboxes/implementVector_parallel_tools.Fortran > `echo implementVector_parallel_tools.Fortran|sed -e "s/implement/complex/" -e "s/Fortran/F90/"`
	@sed -e "s/TAG_/REAL/g" -e "s/tag_/real/g"  < ../toolboxes/implementparallel_tools.Fortran > `echo implementparallel_tools.Fortran|sed -e "s/implement/real/" -e "s/Fortran/F90/"`
	@sed -e "s/TAG_/COMPLEX/g" -e "s/tag_/cmplx/g"  < ../toolboxes/implementparallel_tools.Fortran > `echo implementparallel_tools.Fortran|sed -e "s/implement/complex/" -e "s/Fortran/F90/"`
	@sed -e "s/TAG_/REAL/g" -e "s/tag_/real/g"  < ../layout/implementdatalayout.Fortran > `echo implementdatalayout.Fortran|sed -e "s/implement/real/" -e "s/Fortran/F90/"`
	@sed -e "s/TAG_/COMPLEX/g" -e "s/tag_/cmplx/g"  < ../layout/implementdatalayout.Fortran > `echo implementdatalayout.Fortran|sed -e "s/implement/complex/" -e "s/Fortran/F90/"`
	@sed -e "s/TAG_/REAL/g" -e "s/tag_/real/g"  < ../layout/implementVector_datalayout.Fortran > `echo implementVector_datalayout.Fortran|sed -e "s/implement/real/" -e "s/Fortran/F90/"`
	@sed -e "s/TAG_/COMPLEX/g" -e "s/tag_/cmplx/g"  < ../layout/implementVector_datalayout.Fortran > `echo implementVector_datalayout.Fortran|sed -e "s/implement/complex/" -e "s/Fortran/F90/"`
	@sed -e "s/TAG_/REAL/g" -e "s/tag_/real/g"  < ../input_output/avs_implement.Fortran > `echo avs_implement.Fortran|sed -e "s/implement/real/" -e "s/Fortran/F90/"`
	@sed -e "s/TAG_/COMPLEX/g" -e "s/tag_/cmplx/g"  < ../input_output/avs_implement.Fortran > `echo avs_implement.Fortran|sed -e "s/implement/complex/" -e "s/Fortran/F90/"`
