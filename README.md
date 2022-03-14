# Fiji and Matlab scripts for automated cell tracking and subcapsular sinus segmentation

If you have any questions please contact j.a.pike@bham.ac.uk

* track.groovy performs automated 3D cell tracking using the TrackMate plugin
* produce_ilastik_training_data.ijm creates downsampled training data for training the ilastik pixel classifier
* run_ilastik_clas.ijm runs the ilastik pixel classifier on the full movies
* SCS_postprocessing.m performs postprocessing of cell tracks and SCS segmentations from TrackMate and ilastik respectively

If these scripts are useful please cite our publication

Zhang <em>et al.</em> Recycling of memory B cell supports affinity maturation to antigenic drift (2021) <em>under review</em>.
