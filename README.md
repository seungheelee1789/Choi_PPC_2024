Repository of example data & MATLAB codes (Choi & Lee, 2024)
=============

example data and MATLAB codes of "Locomotion-dependent auditory gating to the parietal cortex mediates flexible multisensory decisions" by Choi and Lee 2024.

Contact to 'ilsong655@gmail.com' if you have any question.

All the analysis were conducted with MATLAB 2019b. Test on MATLAB 2023b also worked.

Add 'Choi Functions' folder with subfolders in the search path.

Please see note at the first part of each MATLAB code for explanation.

For behavior task, data of each trial type are saved in 4-by-4 matrix. See below image for the trial types.
![Trial type](https://github.com/seungheelee1789/Choi_PPC_2024/assets/164326421/b87ed3f0-768e-4605-b457-cc3fe3f341ec)


"Figure 1 and Extended Figure 1" folder
-------------
Example behavior data include the following variables:

EventHz: Sampling rate of locomotion speed measurement.

EventSeq: History of behavior session.
```
Codes in EventSeq
-1st column: vis go (111) / vis no-go (112) / aud go (211) / aud no-go (212) / congruent go (311) / congruent no-go (312) / AgoVnogo (411) / AnogoVgo (412)
  
-2nd column: lick (1) or without lick (2) during the response window
 
-3rd column: hit (1) / miss (2) / correct rejection (3) / false alarm (4)
 
-4th column: correct (1) / incorrect (0)
```
EventSpeed: Locomotion speed (cm/s). -3 to 3s from the stimulus onset. row: trial; column: time point.

LickTime: Animal's licktime.

State: Stationary or moving session

TrialNumber: Number of trial for each trial type.   

<br/>

"PupilData.mat" contains trial-averaged normalized pupil size of 6 mice. Each row contains data of one mouse.
```
Data in the "PupilData.mat"
-1st column: AgoVnogo lick
-2nd column: AgoVnogo without lick
-3rd column: AnogoVgo without lick
-4th column: AnogoVgo lick
```

"Figure 2-3 and Extended Figure 2-4" and "Figure 4 and Extended Figure 5-6"  folder
-------------
Please download the example data here: "[https://drive.google.com/drive/folders/1ajMBy4--RnYhxOwQn_lAy2nMwwr9YR0e?usp=drive_link](https://drive.google.com/drive/folders/1ajMBy4--RnYhxOwQn_lAy2nMwwr9YR0e?usp=sharing)"

Example calcium imaging data include the following variables:

EventF: Calcium activity (dF/F,%). -3 to 3s from the stimulus onset. row: time point; column: neuron; z: trial.

EventHz: Sampling rate of locomotion speed measurement.

EventSpeed: Locomotion speed (cm/s). -3 to 3s from the stimulus onset. row: different trial; column: time point.

ImgHz: Sampling rate of calcium imaging

LocomotionOn(Off)setF: Calcium activity (dF/F,%). -2 to 2s from the locomotion on/offset. row: time point; column: neuron; z: trial.

LocomotionOn(Off)setSpeed: Locomotion speed (cm/s). -2 to 2s from the locomotion on/offset. row: different trial; column: time point.

nTrial: Number of trial for each trialtype.

State: Stationary ('s') or Moving ('s') session

"Figure 5 and Extended Figure 6" folder
-------------
Optogenetic manipulation data include the following variables:

Each cell contains data of single mouse.

Ctrl_LaserTrial: Number of trial for each trial type with laser stimulation during control experiment.

Ctrl_NoLaserTrial: Number of trial for each trial type without laser stimulation during control experiment.

Main_LaserTrial: Number of trial for each trial type with laser stimulation during actual manipulation experiment.

Main_NoLaserTrial: Number of trial for each trial type without laser stimulation during actual manipulation experiment.

"Figure 6" folder
-------------
Example data in AC nVoke experiments include the following variables:

ACXXX_Merged_FData: Calcium activity (dF/F,%). -3 to 3s from the stimulus onset. row: time point; column: neuron; z: trial.
```
Cells in the "ACXXX_Merged_FData"
{1,1} = 5kHz without opto-stimulation; {1,2} = 5kHz with opto-stimulation;
{2,1} = 10kHz without opto-stimulation; {2,2} = 10kHz with opto-stimulation;
```

ImgHz: Sampling rate of calcium imaging

<br/>

Example data in M2AC calcium imaging include the following variables:

Duration: durations of individual locomotion bouts (s).

EventHz: Sampling rate of locomotion speed (Hz)

ImgHz: Sampling rate of calcium imaging (Hz)

LocomotionBoutF: calcium activity (dF/F,%). -2s from the locomotion onset to 2s after the stimulus offset. each cell: one bout; row: time point; column: neuron;

LocomotionOn(Off)setF: Calcium activity (dF/F,%). -2 to 2s from the locomotion on/offset. each cell: one bout; row: time point; column: neuron; z: trial.

LocomotionOn(Off)setSpeed: Locomotion speed (cm/s). -2 to 2s from the locomotion on/offset. each cell: one bout; row: different trial; column: time point.

<br/>

Optogenetic manipulation data (M2AC optogenetic inhibition) include the following variables:

Each cell contains data of single mouse.

Ctrl_LaserTrial: Number of trial for each trial type with laser stimulation during control experiment.

Ctrl_NoLaserTrial: Number of trial for each trial type without laser stimulation during control experiment.

Main_LaserTrial: Number of trial for each trial type with laser stimulation during actual manipulation experiment.

Main_NoLaserTrial: Number of trial for each trial type without laser stimulation during actual manipulation experiment.

