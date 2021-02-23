The folder "IGD+-MaOEA-RS" contains the source code of the IGD+-MaOEA-RS, A-IGD+-MaOEA-RS, R-IGD+-MaOEA-RS, and AR-IGD+-MaOEA-RS.
The folder "Supplementary material" contains the supplementary material of the article "Towards a More Balanced Reference Set Adaptation Method: First Results".

"""
PARAMETERS

H:                  Number of divisions per objective function (Use a positive integer)
flag:               Reference set adaptation method (Use 0, 1, 2, or 3)
problem:            Name of the problem (Use a valid problem name)
m:                  Number of objective functions (Use a positive integer greater than 1)
max_evaluations:    Maximum number of objective function evaluations (Use a non-negative integer)

Supported problems: DTLZ1, DTLZ2, DTLZ3, DTLZ4, DTLZ5, DTLZ6, DTLZ7,
                    DTLZ1_MINUS, DTLZ2_MINUS, DTLZ3_MINUS, DTLZ4_MINUS, DTLZ5_MINUS, DTLZ6_MINUS, DTLZ7_MINUS,
                    IMOP1, IMOP2, IMOP3, IMOP4, IMOP5, IMOP6, IMOP7, and IMOP8.
"""

A computer with Python 3.8 installed is required. 

The program requires the installation of the modules numpy, matplotlib, and scipy. To install the modules using pip write:
	pip3 install numpy
	pip3 install matplotlib
	pip3 install scipy

For running the program, do the following. Go to IGD+-MaOEA-RS/ and write:

IPython console users:
	%run ./main.py 14 0 DTLZ1 3 50000

Windows users:
	python ./main.py 14 0 DTLZ1 3 50000

Linux users:
	python3 ./main.py 14 0 DTLZ1 3 50000

The output files are generated in IGD+-MaOEA-RS/Results/

