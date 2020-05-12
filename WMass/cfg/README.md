## HOW TO RUN COMMANDS:

make sure that all the booleans at the beginning of the run_wmass_cfg.py file
are set according to your needs and dreams and wishes.

test some events to **really** make sure:
```
heppy <name> run_wmass_cfg.py -N 5000 
```

where <name> will be the name of the test output directory with the file and
directory structure
Can also add some options with -o key=value for example -o test=<testname>, where <testname> is defined in the cfg to activate some specific action

### then submit a production to the batch with condor:

```
heppy_batch.py run_wmass_cfg.py -o <localoutputdir> -r /store/cmst3/group/wmass/w-mass-13TeV/ntuples/WHATEVER -b 'run_condor_simple.sh -t 3000 ./batchScript.sh'
```

where <localoutputdir> is a directory on /afs/ which will have the structure and all the .url files (suggest passing the absolute path to it)
and WHATEVER is a useful, clear, and complete description of what the production has saved.
the time passed to -t is in minutes (converted by the script into seconds for condor)
can pass additional options to heppy_batch.py with --option (it works as -o above to pass some key=value pairs to the cfg)

Note:
When running ntuplizer in batch, it might issue some error messages like the following
```
Command (['xrd', 'eoscms.cern.ch', 'existdir', '/eos/cms/store/cmst3/group/wmass/mciprian/heppyNtuples/TREES_ZPILOT_94X_TEST']) failed with return code: -6

creating  root://eoscms.cern.ch//eos/cms/store/cmst3/group/wmass/mciprian/heppyNtuples/TREES_ZPILOT_94X_TEST
Command (['eos', 'ls', '/eos/cms/store/cmst3/group/wmass/mciprian/heppyNtuples/TREES_ZPILOT_94X_TEST']) failed with return code: 2

```
However, if the folders where created and the jobs are running you can ignore them.

When the jobs are done, you have to merge the different Chunks of each process.
To do it, go to WMass/python/postprocessing/condorMerging/ and run the following command

```
python mergeTrees.py -d <localoutputdir> -o /eos/cms/store/<path-to-trees-on-eos>/ --strict --dryRun > trash.txt
```

--strict perform a check on the tree and other utility files.
--dryRun does not run the jobs, allowing you to first check if and which folders were bad (suggest redirecting output to a txt file as in the example)
If some folders were bad, the output will tell you the command to run for resubmission
Otherwise, remove --dryRun and run again (can also remove --strict at that point, and clearly do not redirect the output to a file)
If you are running on MC and only a small fraction of trees is missing, you can decide to ignore them: in this case, before running without --dryRun, move the bad folders from <localoutputdir> into another location (better not to delete them, you may want to recover them later).

Once this step is complete, you will find N additional folders on eos in the same location were the initial ones where, named as the original one with a postfix _partXX,
Note that the XX will go from 1 to N (usually a different integer for each different process), but it doesn't have a relation with the actual MC files as seen in DAS. 
Therefore, in case you are making ntuples from the same sample but in two steps (e.g. processing half the files in DAS on one day and the rest in another day), it is suggested producing the ntuples in two separate folders, and put them together only at the end (you'll have to rename some parts, as, part_XX will always start from 1 for each production)

Next step is to produce the friend trees to add additional variables (there is a readme in WMass/python/plotter/ to explain how to do it)
Moreover, once the friends are done, you may also want to run a skim on these ntuples to run faster on them
