STEPS FOR BASEBALL ANALYSIS
 0  Download latest version of lahman database.
 1  Preprocess initial data
 2  Run MCMC sampler (run 3-5 consecutive chains saving as output1, output2, etc.)
 3  Preprocess prediction data
 4  Predict homeruns  (run 3-5 consecutive chains saving as output1, output2, etc.)
 5  Compile predictions and make trace plots.

Notes:
- In step 2, you only need to change the save() filename on the very last line.
- In step 4, you need to change a load() and a save() command.  The lines are marked by a #*#.
- When I updated the results for a talk a few years ago (i.e., same model but on more recent data), I had to insert a kluge for the Beta draws in Step 2 (it is clearly marked).  In particular, for the conditional MLE proposal was sometimes off for RF bc all right fielders Elite indicators were either 0 or 1 on a given draw thereby causing problems fitting the MLE.  For that dataset, there were only 44 RF player-seasons (most were listed as OF) but this was not a problem in our paper.  The kluge lines are clearly marked and can be removed.
