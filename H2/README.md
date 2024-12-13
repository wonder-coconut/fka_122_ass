# Homework 2 part B - FKA122 Computational Physics
### Compilation
You can compile the script using the makefile, just use the `make` command. To remove all objects, use the `make clean` command. The program is compiled into the `H2` object

### Execution
The script, as of the moment of updating this document, contains two functions which are sub-problems of task 1 (sampling and correlation). Each of these tasks can be executed using the first command line argument.The resultant execution format varies based on the task.


#### Task 1, subtask 1
The sampling subtask uses MCMC sampling to sample the trial wavefunction for the He atom system. Walker updates are done defined as:
$new\_position = original\_position + delta * (uniform\_rv - 0.5)$. The delta value chosen for this update such that the acceptance ratio comes close to 0.5 is 0.972.

Execution format:
```
./H2 1 <number of samples> <delta value>
```
An ideal example execution would be `./H2 1 1000000 0.972`

#### Task 1, subtask 2
The correlation subtask requires no further arguments.

Execution format:
```
./H2 2
```


#### Task 2, subtask 1
Similar to T1S1, this is a sampling subroutine, with the only difference of having an "unlikely" initial starting point. The syntax of execution is identical to T1S1.
```
./H2 3 <n samples> <delta value>
```

#### Task 2, subtask 2 & 3
This part handles the statistical inefficiency calculations, and similar to T1S2 requires no parameters.

Autocorrelation method:
```
./H2 4
```
Block averaging method:
```
./H2 5
```

#### Task 3
T3 takes the same arguments as T1, except here, the script runs for multiple alpha values

```
./H2 6 <n samples> <delta value>
```

#### Task 4
T4 takes additional arguments for the maximum p value, the initial alpha value, and the beta parameter.

```
./H2 7 <n samples> <delta value> <p_limit> <alpha_init> <beta>
```

#### Task 5
T5 is basically identical T1, and so are its parameters.
```
./H2 8 <n samples> <delta value>
```
The remainder of the project involves visualization by python, which I assume to be self-evident in their execution.