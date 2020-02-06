# Figures for point-to-point task #

This collection is associated with the article "De novo learning
versus adaptation of continuous control in a manual tracking task."
This code generates all figures related to the sum-of-sinusoids
tracking task. Raw data, organized by subject number, is contained in
the folder "Data." Within each subject's data folder is folders
corresponding to different blocks of the experiment:

    baseline: reaches without rotation or mirror-reversal
    pert1: first point-to-point block under either perturbation
    pert2: second block under either perturbation
    pert3: third block under either perturbation

Data from individual trials are stored in separate files within each
folder for a block. Also included is the "tFile" which records the
target position on each trial. The columns indicate the following:

    column 1: trial number
    column 2: target x-position
    column 3: target y-position
    columns 4-5: extraneous information that is not used for analysis

All analyses can be performed and all figures generated by simply
running loadBimanualSkillData.m. All other .m files are functions that
are used by main.m for analysis or plotting. Briefly, these functions
do the following:

    loadSubjData.m: extract raw data from data files
    plot_direction.m: plots reach-direction error
    processData.m: analyzes raw data
    savgolayFilt.m: Savitzky-Golay filter to smooth trajectories