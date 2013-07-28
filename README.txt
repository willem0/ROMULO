README.txt

  This MATLAB code is released under the Gnu Public License (GPL). For more information, 
  see : http://www.gnu.org/copyleft/gpl.html

Please direct all questions to "William Mantzel" <willem@gatech.edu>.

The main scripts of interest are:
--- mscmfp_batch.m:       the batch script which configures the simulation and runs...
--- mscmfp_simulate.m:    the script that actually runs the thousand simulations-per-line, saving the results in the data/ directory
--- mscmfp_figGen.m:      generates the output figures from the saved data files

The supporting functions are:
- main_modes.m:           generate the modes of the Pekeris waveguide, saving the result in the states/ folder
- normal_modes.m:         helper function for main_modes.m
- greens_mode.m:          generate the Green's function (complex frequency response) between a source location and set of receivers
- greensG_mode.m:         generate the Green's function between a series of potential source locations (e.g., grid) and a set of receivers
- mscmfp_global.m:        simulate of the global (as opposed to greedy) least-squares localization of a pair of sources
- mscmfp_clean.m:         simulate of the CLEAN localization method of Song et al.
- mscmfp_lipchitz.m:      compute the upper and lower Lipchitz bounds for the Green's function with respect to source location
- minimax_assignment.m:   choose the source estimate relabeling that minimizes the maximum distance to the actual sources
- munkres.m:              helper function for minimax_assignment.m that utilizes the Hungarian algorithm for assignment
- mscmfp_sigmoid.m:       utility function that computes a sigmoid function
- mscmfp_gsmooth.m:       utility function that performs Gaussian smoothing
- normByCols.m:           utility function that computes the norm of each column in a matrix
