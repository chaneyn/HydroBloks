import os
#Compile NOAH
os.chdir('src/Utility_routines')
os.system('make clean')
os.system('make')
os.chdir('../Noah')
os.system('make clean')
os.system('make')
os.system('cat pyNoahMP.F > pyNoahMP.f90')
#Create signature file
#cmd = 'f2py pyNoahMP.f90 -h NoahMP.pyf -m NoahMP --overwrite-signature'
#os.system(cmd)
#Define subroutine to import 
subroutines = 'initialize'
#Create driver library
#cmd = 'f2py -c NoahMP.pyf'
cmd = 'f2py -c pyNoahMP.f90 *.o ../Utility_routines/*.o -m NoahMP -I../Utility_routines -L/lib64 -L/usr/lib64 --f90flags="-w -fopenmp -g -Werror -fmodule-private -fimplicit-none -fbounds-check -fcheck=array-temps,bounds,do,mem,pointer" only: initialize update :'
#cmd = 'f2py --debug -lgomp -c NoahMP.pyf *.o ../Utility_routines/*.o --fcompiler=gnu95 -L/lib64 -L/usr/lib64 --f90flags="-w -fopenmp -g -Werror -fmodule-private -fimplicit-none -fbounds-check -fcheck=array-temps,bounds,do,mem,pointer"'
#cmd = 'f2py -c only: %s --debug -lgomp -c NoahMP.pyf *.o ../Utility_routines/*.o --fcompiler=gnu95 -L/lib64 -L/usr/lib64 --f90flags="-w -fopenmp -g -Werror -fmodule-private -fimplicit-none -fbounds-check -fcheck=array-temps,bounds,do,mem,pointer"' % subroutines
#cmd = 'f2py -m NoahMP --debug -lgomp -c only: run_model run_model_cell : pyNoahMP.f90 *.o --fcompiler=gnu95 -I../Utility_routines -L/lib64 -L/usr/lib64 --f90flags="-w -fopenmp -g -Werror -fmodule-private -fimplicit-none -fbounds-check -fcheck=array-temps,bounds,do,mem,pointer -ffree-form  -ffree-line-length-none"'
#cmd = 'f2py -m NoahMP --debug -lgomp -c pyNoahMP.f90 *.o ../Utility_routines/*.o --fcompiler=gnu95 -I../Utility_routines -L/lib64 -L/usr/lib64 --f90flags="-w -fopenmp -g -Werror -fmodule-private -fimplicit-none -fbounds-check -fcheck=array-temps,bounds,do,mem,pointer -ffree-form  -ffree-line-length-none"'
#print(cmd)
#cmd = 'f2py -m NoahMP --debug -lgfortran -lgomp -c NoahMP.pyf *.o ../Utility_routines/*.o only: initialize_parameters run_model : --fcompiler=gnu95 -L/lib64 -L/usr/lib64 --f90flags="-w -fopenmp -g -Werror -fmodule-private -fimplicit-none -fbounds-check -fcheck=array-temps,bounds,do,mem,pointer -ffree-form  -ffree-line-length-none"'
#cmd = 'f2py -m NoahMP --debug -lgomp -c NoahMP.pyf *.o ../Utility_routines/*.o only: initialize_parameters run_model : --fcompiler=gnu95 -L/lib64 -L/usr/lib64 --f90flags="-w -fopenmp -g -Werror -fmodule-private -fimplicit-none -fbounds-check -fcheck=array-temps,bounds,do,mem,pointer"'
#exit()
#cmd = 'f2py --debug -c NoahMP.pyf *.o ../Utility_routines/*.o --fcompiler=gnu95 -L/lib64 -L/usr/lib64 --f90flags="-w -g -Werror -fmodule-private -fimplicit-none -fbounds-check -fcheck=array-temps,bounds,do,mem,pointer"'

os.system(cmd)
#Clean up both directories
#os.system('mv NoahMP.*.so ../../NoahMP.so')
os.system('mv NoahMP*so ../../NoahMP.so')
