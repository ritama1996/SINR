all :
	gfortran -c read_file_mod.f03
	gfortran -c set_parameter_mod.f03
	gfortran -c set_input_mod.f03
	gfortran -c get_force_mod.f03
	gfortran -c update_variables_mod.f03
	gfortran -c mts_method_mod.f03 
	gfortran -c molecular_dynamics_mod.f03
	gfortran -c sinr_algorithm.f03
	gfortran -o sinr.exe -g read_file_mod.o set_parameter_mod.o set_input_mod.o get_force_mod.o update_variables_mod.o mts_method_mod.o molecular_dynamics_mod.o sinr_algorithm.o
