temporary fix:
run this before executing ./H1

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/wonder_coconut/.local/lib/
export LD_LIBRARY_PATH

makefile commands:

    make : compile all objects

    make clean : delete all compiled objects

execution:

    "H1" is the main program object
    H1 takes three command line arguments, depending on the task
    Arg 1 is the choice for the task
    Arg 2 can either be:
        1. timestep, for choices 2-5
        2. n_max, for static factor calculation (choic 6)
    Arg 3 is the t_max for choices 2-5