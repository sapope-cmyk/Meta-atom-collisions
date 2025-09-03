# Meta-atom-collisions
These are the Matlab files to create the data for the surface and time series plots for the journal paper "Meta-atom based two-sphere Newton's cradle" https://doi.org/10.1103/2whk-jc4b

To generate the data for the surface plots:
1. Run Model_v3.m to setup the outer sphere poaranmeters and environment
2. Run Code_to_generate_data_for_surface_plots.m to generate the surface plot data (output variables listed in the file)

To simulate the Meta-atom based two-sphere Newton's cradle and generate the time series data:
1. Run Model_v3.m to setup the outer sphere poaranmeters and environment
2. Run Code_to_set_initial_conditions_for_the_time_series_simulation.m (set the required resonator frequency and initial conditions in this file - values are contained for the simulations presented in the paper - comment out as necessary)
3. Run Meta_atom_two_sphere_Newtons_cradle_simulation.slx to start the simualtion and generate the time series data
