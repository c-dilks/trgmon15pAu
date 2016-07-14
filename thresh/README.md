Summary
-------
Execute `getThresh` to run the following commands:
- `root -b -q RunIndex.C`
- `loopThresh`
- `root -b -q BuildThreshTree.C`
- `root -b -q DrawThreshPlots.C`
- `copyToProtected`
 

Procedure
---------
These scripts pull trigger threshold information from trigger version webpage
and builds a table of run numbers together with corresponding trigger 
thresholds; this table is named `thresh.dat`

1. Obtain `runlog.root` file from trgmon results

2. Execute `RunIndex.C` to generate a run index file, `index.dat`
   - The columns are `[ run# | index# | sqrt(s) | download boolean ]`
   - To do this, run `RunIndex.C("runlog.root",sqrt(s))`
   - `doubleRunIndex` will do this for two separate root files,
     `runlog200.root` and `runlog500.root` (from `trgmon12`) and
     concatenate an `index.dat` file
   - download boolean is `true` if info needs to be downloaded

3. Execute `loopThresh`
 - This will build `thresh.dat`, a file containing trigger thresholds
   for the run numbers listed in `index.dat`
 - This runs `triggerThresholds` and `appendThreshFile` (see subroutines below)
   for all runs listed in `index.dat`; in principle, to update plots, one
   can generate a new `index.dat` file containing new runs and compare it to
   the old one so that only new threshold information is downloaded; then
   a new `thresh.root` file can be created and hadded to the old one
 - temporary `trg_*` and `th_*` files used in `triggerThresholds` will be then
   relocated to subdirectories, but are probably not useful after building
   `thresh.dat` and can be deleted 
 - Columns of `thresh.dat` (left to right) are given in `appendThreshFile`
    

4. Execute `BuildThreshTree.C` to build a tree from `thresh.dat`, 
   which will be written to `thresh.root`

5. Execute `DrawThreshPlots.C` (use `-b -q`) to draw threshold plots vs. run
   index; these data are written to `th_plots.pdf`

6. Execute `copyToProtected` to update the link from the webpage



SUBROUTINES
-----------
To obtain trigger thresholds for a single run, execute
  `triggerThresholds [run#]`
   - two files are created; `trg_[run#]` is used to obtain triggerid
     and trigger version numbers which are then subsequently used
     to retrieve the appropriate trigger thresholds, written to
     `th_[run#]`

To append the relevant trigger numbers to `thresh.dat`, execute
  `appendThreshFile [run#] [index#] [sqrt(s)]`
