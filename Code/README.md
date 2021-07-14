# GeometryGlycoproteinEvolutionSimulations


The following files were used to create Lammps simulations that compute the mean first passage time of the antibody to one of the epitopes on the surface of hemagglutinin.

“Virus_Flu_Spikes${ SpikeN}_Cond${Cond}”: Structure file and input to Lammps. Contains the coordinates of the beads in the system. This is the virus model presenting “SpikeN” HA molecules on its surface. Each of the 6 conditions (“Cond”) corresponds to a different initial configuration/position of the antibody molecule.
To submit the simulations use the script “submit_multiple_jobs.sh”. This script is designed to submit multiple simulations with for different conditions.

“Virus_Corona_Spikes${SpikeN}_Cond${Cond}”: Structure file and input to Lammps. Contains the coordinates of the beads in the system. This is the virus model presenting “SpikeN” S proteins on its surface. Each of the 4 conditions (“Cond”) corresponds to a different initial configuration/position of the antibody molecule.
To submit the simulations use the script “submit_multiple_jobs.sh”. This script is designed to submit multiple simulations with for different conditions.

“SpikeNum”: the number of spikes (HA molecules or S proteins) on the surface of the virus.

“cond” initial position of the antibody.

“idumNum”: seed number for the random number generator of Lammps.

“epitope”: the target residue (epitope) on the surface of HA. There is a total of 228 epitopes.

“lammps_flu.pbs”: setup file for the flu virus simulation submission to a slurm cluster.

“FluSim_epitope.in”, “CoronaSim_epitope.in”: setup file for the Lammps flu and corona virus simulations. They contains all the simulation’s variables.

The simulation outputs a file “SpikeSim_idum_${idum}_virushaSpikes${SpikeN}_Bond_Ep_${epitopenum}_Cond${Cond}”
When the variables here correspond to the input variables. This is a text file detailing whether or not the antibody molecule arms (Fab) touch their respective epitope. From this file I extract the on-rate for encounter (inverse of the mean first passage time).

“ReadBondInfo_MFPT_Epitope_RatesAllClones.m”: Matlab script to analyze the output files “SpikeSim_idum_${idum}_virushaSpikes${SpikeN}_Bond_Ep_${epitopenum}_Cond${Cond}” and extract the on-rate of the first and second antibody arms to the epitope.

The second output of the simulation is a file of the format:
SpikeSim_idum_${idum}virushaSpikes${SpikeN}_Ep_${epitopenum}_Cond${Cond}

It allows visualizing the simulation with a software such as Ovito

Attached is an example “Virus_Flu_Spikes56_Ep_158_Cond6”


The code is related to the manuscript [**Viral surface geometry shapes influenza and coronavirus spike evolution**](https://www.biorxiv.org/content/10.1101/2020.10.20.347641v1)

To look at the analysis to the time-dependent evolution of SARS-CoV-2 go [**here**](https://amitaiassaf.github.io/SpikeGeometry/SARSCoV2EvoT.html)