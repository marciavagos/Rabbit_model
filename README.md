This model is written in C++. To compile the source code, use eg. the following command:

    g++ main.cc RAMM.cc -O3 -o RAMM.out

RAMM.cc contains the code of the model itself, with the ion current formations, balance equations, etc.
main.cc contains helper code to define user-specified input parameters and stdout output.

To run a simulation in the bash shell, run the following command:

    ./RAMM.out Parameters.txt output CELL

where <output> is to be replaced with the desired output file name. 
The file Parameters.txt is a settings file that contains the user defined parameters, 
and has to be specified. In this file you can define simulation settings such as pacing cyle length,
duration, stimulation current, voltage clamp protocol, etc.

An examplary parameter file is provided in Parameters_baseline.txt, which runs a simulation
with the baseline (publish) model parameters.

To run a population of models, you can execute the bash script provided in the population
folder (<population.sh>), which makes use of the <Parameters_batch.txt> file to run simulations
with different maximum conductances scalings (relative to baseline parameters). 

The file <population/params_CTL.txt> rovides 17 different maximum conductance scalings,
the first line corresponding to the baseline values, and the other lines corresponding to the 16 models published in the [Vagos et al., 2020](https://www.frontiersin.org/articles/10.3389/fphys.2020.556156/full) paper.

