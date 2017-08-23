---
title: "smt change log"
output: html_document
---

## 0.4.0 - 2017-08-21
### Added
- Created a kernel density masking algorithm that uses a kernel density calculation along with masking application algorithms to create a trackll mask based on density
- Added features such as automatically removing edge clusters, separating clusters, and model building/loading into the .densityMaskTracks and densityMaskTracks
- Created plotPoints() and plotLines() and their trackl counterparts to simply plot all trajectories
- Made .readDiaSessions compatible with newer versions of Diatrack
- Made compareFolder compatible with createTrackll and all recognized trackll input file types

### Deprecated

### Removed
- Masking feature base on image masks removed and made into a new function maskTracks
- Merging feature removed and made into a new function mergeTracks

## 0.3.9 - 2017-07-22
### Added
- Created a .readDiaSessions algorithm to process a single Diatrack session file into a track list with the ability to calculate absolute coordinates and record frames.
- Implemented .readDiaSessions into readDiaSessions to process multiple session files into a list of track lists in parallel with the ability to merge, calculate absolute coordinates, mask, and record frames.
- Created a .readSlimFast algorithm to process a single SlimFast file into a track list with the ability to calculate absolute coordinates and record frames.
- Implemented .readSimFast into readSlimFast to process multiple session files into a list of track lists in parallel with the ability to merge, calculate absolute coordinates, mask, and record frames.
- Added the feature in .readDiatrack and readDiatrack to artificially generate frame records in consecutive order after the recorded start frame.
- Added the feature in .readParticleTracker and readParticleTracker to collect frame records from the input file.
- Created a createTrackll script that calls the different read functions to create a list of track lists from any desired input with the option to interact with a selection memu. All read parameters have been preserved and the ability to add more read fucntions in the future has been implemented.
- Created a .exportRowWise script to export any track list (linked, processed, or from any originating source) into a ImageJ/MOSAIC style row-wise .csv file output.
- Implemented .exportRowWise into a exportTrackll script to export a list of track lists into multiple .csv save files. Added the ability to optimize multiple cores in parallel. Exporting a single file or multiple files has the ability to be read back in with the .readParticleTracker/readParticleTracker/createTrackll scripts.
- Created a .linkSkippedFrames script to link any seemingly skipped frames in a track list given a tolerance level in pixels and the maximum number of frames it could skip. Adds the number of skips to the track names.
- Implemented .linkSkippedFrames into a linkSkippedFrames script to link track lists in a list of track lists. Added the ability to optimize multiple cores in parallel.
- Added a getStartFrame helper function to simply return a single track's start frame using its name and "." delimiters.
- Added a getTrackFileName helper functio to return a track's abbreviated five character subname using its name and "." delimiters.
- Added a removeFrameRecord helper function to remove the frame record from the fourth column of any track list in order to preserve backwards compatibility with existing analysis functions.
- Documented all added features and changes into their respective helpers using Roxygen/devtools.

### Changed
- Replaced the computation of absolute coordinates of a track into the helper function abTrack and implemented in all reading algorithms.
- Changed the way .readParticleTracker reads in the file subname to ignore any ".tif" before the .csv extention and take the last five characters of the file name before either of these extensions.
- All read functions adjusted to have the same parameters of calcualting the absolute coordinates and recording the frames. The ability to merge, mask, and optimize cores were kept uniform as well.
- Made the standard print output of every function uniform in terms of 
spelling, displaying processing and confirmation text, and spacing.

### Deprecated
- Plan on removing code duplication of parallelization and replace it with a parallel function.

### Removed
- Added a .exportColWise script to export track lists into a pseudo-Diatrack style output, but is removed for now (can be found in exportTrackll).
- Tested a findMaxLinks script to find the optimum tolerance and maximum skip number, but is removed for now (can be found in linkSkippedFrames).
