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

The file <population/params_pop.txt> provides 3000 different combinations of 13 selected model parameters
(including  maximum conductance scaling factors) that used in the population of models study described in 
[Vagos et al., 2020](https://www.frontiersin.org/articles/10.3389/fphys.2020.556156/full) paper.
Of these, 16 models where selected that represent physiological data from the rabbit atrial cardiomyocyte 
(more details in the paper). Those 16 models correspond to the indices in the <array> variable in the <population/population.sh> script.
The first line corresponding to the baseline values, and the other lines corresponding to the 16 models 
published in the paper.
The file <population/params_CTL.txt> contain the parameter sets of baseline + 16 models only.

To run a batch population of models, you can execute the bash script provided (<population/population.sh>),
which makes use of the  <Parameters_batch.txt> file to run simulations with different parameter settings.
You can either re-run the whole population of 3000 models (eg., to  experiment with other model parameter changes),
or just the 16 selected model by specifying the correct params file (<params_pop.txt> or <params_CTL.txt>) and the 
corresponding indices in the <array> variable.
 
Additional simulated data of this  population of 16 models can be found [here](https://github.com/marciavagos/rabbit_model_datasets).

