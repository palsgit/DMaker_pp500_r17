# Dmaker_pp500

Compile on STAR's RCAF:
```sh
starver SL22b
cons
```

Run D0 analysis locally on few picoDsts with:
```sh
root -l -b -q StRoot/macros/loadSharedHFLibraries.C StRoot/macros/runPicoD0AnaMakerLocal.C
```
Mixed-event background pairs of kaons and pions:
```sh 
root -l -b -q StRoot/macros/loadSharedHFLibraries.C StRoot/macros/runPicoMixedEventLocal.C
```
QA of picoDst data:
```sh
root -l -b -q StRoot/macros/loadSharedHFLibraries.C StRoot/macros/runQAAnaMakerLocal.C
```


Actual file list is loaded from:
```sh
picoLists/runs_path_all.list
```
Actual file list for **local tests** is loaded from:
```sh
picoLists/runs_local_test.list
```
After creating your directories and setting your name in run*.sh scripts, run analysis on farm by run-scripts:
```sh
runJob.sh
runJobMixed.sh
runJobQA.sh
```
Merging og job output files with script
```sh
./merge_new.sh production/
```
Analysis of job outputs in folder /scripts
```sh
Danalysis_pp.s
```

# DMaker_pp500
