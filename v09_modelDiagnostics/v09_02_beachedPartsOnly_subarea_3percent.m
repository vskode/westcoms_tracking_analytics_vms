% use this script to convert the final positions of x and y to only beached coordinates

clear 
close all 

% load model output 
load ../cluster_output/01_ModelRuns_ClydePaper/v09_02_3percent/v09_02_output_edited_3percent.mat  

% count how many ashore and inwater at the end of the model run 

onbeach = numel(ashore(ashore));
inwater = numel(ashore(~ashore));

% concatenate the positions and logical variable of ashore 
fate = [ashore, x, y]; 

% select only points that equal 1 (beached)
onlybeached = fate(fate == 1,:); 

% use only x and y 
onlybeached = onlybeached(:,[2 3]); 

% save each only individually 

beachedx = onlybeached(:,1);
beachedy = onlybeached(:,2); 












