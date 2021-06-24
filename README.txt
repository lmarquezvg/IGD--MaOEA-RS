Documentation of the module main.py used for the experimentation of: Towards a 
More Balanced Reference Set Adaptation Method: First Results.

NAME
       main.py - test a multi-objective evolutionary algorithm (MOEA)

SYNOPSIS
       main.py H FLAG MOP OBJS EVALS
       main.py OPTION

DESCRIPTION
       This module is used to test different versions of the IGD+-based MOEA 
       with uniform reference set on a selected multi-objective problem (MOP) 
       with a given number of objective functions. The required arguments are 
       described as follows:

       H
              It must be an integer greater than zero. It represents the number 
              of divisions per objective function used for the generation of a 
              uniform reference set. The values 99 and 14 are selected for MOPs 
              with 2 and 3 objective functions, respectively.

       FLAG
              It is used to determine which reference set adaptation method 
              will be employed during the execution of the IGD+-based MOEA with 
              uniform reference set. The value zero means that no adaptation 
              method is used, and the values one, two, and three correspond to 
              the adaptation methods from A-NSGA-III, RVEA*, and AR-MOEA, 
              respectively.

       MOP
              It must be a valid MOP name. The valid MOP names are: DTLZ1, 
              DTLZ2, DTLZ3, DTLZ4, DTLZ5, DTLZ6, DTLZ7, DTLZ1_MINUS, 
              DTLZ2_MINUS, DTLZ3_MINUS, DTLZ4_MINUS, DTLZ5_MINUS, DTLZ6_MINUS, 
              DTLZ7_MINUS, IMOP1, IMOP2, IMOP3, IMOP4, IMOP5, IMOP6, IMOP7, and 
              IMOP8.

       OBJS
              It must be an integer greater than one. It represents the number 
              of objective functions of the MOP. A value of two is used for 
              IMOP1-IMOP3. While a value of three is used for DTLZ1-DTLZ7, 
              DTLZ1_MINUS-DTLZ7_MINUS, and IMOP4-IMOP8. 

       EVALS
              It must be an integer greater than or equal to zero. It 
              represents the maximum number of objective function evaluations 
              for the execution of the MOEA. In our experiments, a value of 
              50000 is used in all cases. 

       The following option can be used:

       --help 
              Display this help and exit.

REQUIREMENTS
       A computer with the installation of Python 3.8 is needed. The modules 
       numpy, matplotlib, and scipy are required.

EXAMPLE
       For running the module main.py, go to IGD+-MaOEA-RS/ and write:

       IPython console users:
              %run main.py 14 0 DTLZ1 3 50000

       Windows users:
              python main.py 14 0 DTLZ1 3 50000

       Linux users:
              python3 main.py 14 0 DTLZ1 3 50000

       The previous line executes the module main.py to test the IGD+-based 
       MOEA with uniform reference set using no adaptation method on the DTLZ1 
       with three objective functions. A uniform reference set with fourteen 
       divisions per objective function is used. The maximum number of 
       objective function evaluations is set to 50000.

RESULTS
       On success, the output files containing the approximation set are 
       generated in IGD+-MaOEA-RS/Results/.

CONTENTS
       The folder "IGD+-MaOEA-RS" contains the source code of the module 
       main.py. The folder "Supplementary material" contains the supplementary 
       material of the paper: Towards a More Balanced Reference Set Adaptation 
       Method: First Results.

