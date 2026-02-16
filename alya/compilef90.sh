#mpif90 -cpp -DUSEMPIF08 unitt_alya_with_another_code_enhanced_no_time_loop.f90 -o Alya.x
#mpif90 -cpp -DUSEMPIF08 alya_time_loop.f90 -o Alya.x
#mpif90 -cpp -DUSEMPIF08 alya_all2all.f90 -o Alya.x
mpif90 -cpp -DUSEMPIF08 alya_all2all_time_loop.f90 -o Alya.x
