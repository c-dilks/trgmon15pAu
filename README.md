This is the code that builds `www.star.bnl.gov/protected/spin/dilks/trgmon15pAu/200/trigger.html`

Procedure
---------
1. Execute `getTrgEv` to download all needed data from the `RunLog` database
   - `trigger.dat` is the main file that will be used
   - there are if statements in `getTrgEv` which select what goes into 
     `trigger.dat`, such as requirents for nonzero no. of `BBC` events and 
     `JP1` events
   - an index is assigned to each run number, starting from the first one
     that is appended to `trigger.dat`; this index is an integer assigned
     to each run number so plots of things vs. run number are plotted 
     vs. the index instead so they are easier to read
2. Execute `mkHTML` to build the HTML pages
   - this runs `trgMon.C`, which is important and described below   
   - each page is written using other `csh` scripts
   - data files and HTML files will be copied to `./200` or `./500`
     depending on `sqrt(s)` being used
   - if `RunLog` SQL database is inaccessible, all the `dat` files that
     are used to build the HTML page are always maintained in
     `./200` or `./500`; check the SQL server port map to see if
     the database has been moved and update `getTrgEv` accordingly
3. Execute `copyToProtected`, which will copy all the files
   to the proper directory for HTML access


`trgMon.C`
----------
This script attempts to determine which runs are good to keep, 
by comparing the trigger rate of a run to the average trigger rate
within an epoch. An epoch is a duration of time where the trigger
rate is approximately constant; the epoch separators are inserted
manually in `trgMon.C` for each trigger. 

Trigger rate is computed as: `(#events * prescale)/(#min. bias * min. bias prescale)`
 - The trigger rates are all drawn in `TCanvas *ev_canv`

Then the trigger rate is normalised by the epoch mean in each epoch;
 - These normalised plots are drawn in `TCanvas *nv_canv`

Then cutoffs, defined under the epoch cuts, are used to decide whether a run
should be marked as good or bad; the lower cutoff is called `qa_cut_low` and is less than
unity and the higher cutoff is called `qa_cut_high` and is greater than unity
 - The cutoff lines define a region around unity in the epoch mean normalised 
   trigger rate plots (cno plots) for which runs are marked as good

A run is only marked as good if the duration is greater than 3 minutes (`time_cut`)

Sometimes, a run will cause anomolously high or low trigger rates and radically
throw off the epoch mean; for these runs, we can manually mark them as bad
and remove them from the computation of the epoch mean; this is done via the `exclusion_list`

Any other bad runs to be manually excluded are added to the variable `exclusion_list`

Trigger ratios are plotted to `TCanvas *rat_canv`

Luminosity plots are drawn automatically, to `TCanvas *lum_canv`

Prescale vs. run index are drawn to `TCanvas *prescale_canv`

Everything is written to `runlog.root`
