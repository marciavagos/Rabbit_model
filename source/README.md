This model is written in C++. To compile the source code, use eg. the following command:

    g++ main.cc RAMM.cc -O3 -o RAMM.out

To run a simulation in the bash shell, run the following command:

    ./RAMM.out Parameters.txt output CELL

`RAMM.cc` contains the code of the model itself, with the ion current formations, balance equations, etc. 

`main.cc` contains helper code to define user-specified input parameters and stdout output.
