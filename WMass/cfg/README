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

where <localoutputdir> is a directory which will have the structure and all the .url files
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
