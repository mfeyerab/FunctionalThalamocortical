
TrainingSST = readtable(...
    "D:\Cloud\OneDrive - The University of Western Ontario\Rebutal\Reanalysis\SST\SST_Classification.xlsx");
IC = findgroups(TrainingSST.IC);
MType = findgroups(TrainingSST.Var5);

%%

[p,h, stats] = ranksum(TrainingSST.Tau_fit_(MType==1 & IC==1), ...
    TrainingSST.Tau_fit_(MType==2 & IC==1))
[p,h, stats] = ranksum(TrainingSST.Rin_hd_(MType==1 & IC==1), ...
    TrainingSST.Rin_hd_(MType==2 & IC==1))
%% Classification Cesium
input = [TrainingSST{IC==1,[7,6]}];
label = MType(IC==1);
UnknownSST = readtable("D:\Cloud\OneDrive - The University of Western Ontario\Rebutal\Reanalysis\SST\SST_prediction.xlsx");
reps =100;

threshold = 0.65;
PredMat = zeros(height(UnknownSST),reps);
Perf = zeros(1,reps);

for r=1:reps

Mdl = fitcknn(...
    input, ...
    label, ...
    'Distance', 'Euclidean', ...
    'Exponent', [], ...
    'NumNeighbors', 3, ...
    'DistanceWeight', 'Equal', ...
    'Standardize', true);

partMdl = crossval(Mdl, 'kfold',23);
Perf(1,r) = 1 - kfoldLoss(partMdl, 'LossFun', 'ClassifError');
PredMat(:,r) = predict(Mdl,UnknownSST{:,[5,4]});

end

boxchart(Perf)

%% Getting Majority Label from all repetitions

for r=1:height(UnknownSST)
   temp= tabulate(PredMat(r,:));
    if temp(1,3)==100 
        UnknownSST.newSSTlabels(r) = temp(1,1);
    elseif  temp(1,3)==0
        UnknownSST.newSSTlabels(r) = temp(2,1);
    end
end

%% Making a horizontal bar plot indicating distributions of putative M-types across layers

L1MC = sum(UnknownSST.newSSTlabels(string(UnknownSST.Var3)=="L1")==1);
L1NMC = sum(UnknownSST.newSSTlabels(string(UnknownSST.Var3)=="L1")==2);
L2_3MC = sum(UnknownSST.newSSTlabels(string(UnknownSST.Var3)=="L2/3")==1);
L2_3NMC = sum(UnknownSST.newSSTlabels(string(UnknownSST.Var3)=="L2/3")==2);
L4MC = sum(UnknownSST.newSSTlabels(string(UnknownSST.Var3)=="L4")==1);
L4NMC = sum(UnknownSST.newSSTlabels(string(UnknownSST.Var3)=="L4")==2);
L5AMC = sum(UnknownSST.newSSTlabels(string(UnknownSST.Var3)=="L5A")==1);
L5ANMC = sum(UnknownSST.newSSTlabels(string(UnknownSST.Var3)=="L5A")==2);
L5BMC = sum(UnknownSST.newSSTlabels(string(UnknownSST.Var3)=="L5B")==1);
L5BNMC = sum(UnknownSST.newSSTlabels(string(UnknownSST.Var3)=="L5B")==2);
L6MC = sum(UnknownSST.newSSTlabels(string(UnknownSST.Var3)=="L6")==1);
L6NMC = sum(UnknownSST.newSSTlabels(string(UnknownSST.Var3)=="L6")==2);

barh({'L1', 'L2/3', 'L4', 'L5A', 'L5B', 'L6'}, [...
    L1MC,L1NMC; L2_3MC,L2_3NMC;  L4MC, L4NMC; L5AMC, L5ANMC; ...
    L5BMC, L5BNMC ; L6MC, L6NMC],'stacked')
set(gca, 'YDir', 'reverse');

%% Getting effects of M-type and thalamic nucleus on post synaptic features 

TrainingSST.ID(contains(TrainingSST.ID,'180109T1_S2C1')) = {'180109_S2C1'};% Correction of mismatching labels
Response = readtable("D:\Cloud\OneDrive - The University of Western Ontario\Rebutal\Reanalysis\SST\SST_responses.xlsx");
Response = addvars(Response,nan(height(Response),1));

for r=1:height(Response)
    if any(ismember(string(UnknownSST.Var1),string(Response.Var1{r})))
      Response.Var15(r) = ...
       UnknownSST.newSSTlabels(ismember(string(UnknownSST.Var1),string(Response.Var1{r})));
    elseif any(ismember(string(TrainingSST.ID),string(Response.Var1{r})))
      if  string(TrainingSST.Var5{ismember(string(TrainingSST.ID),string(Response.Var1{r}))}) == "nMC" 
        Response.Var15(r) = 2;
      else
        Response.Var15(r) = 1;
      end          
    else
        disp(Response.Var1{r})
    end
end

Response.Var15(isnan(Response.Var15)) = 1; % Pooling L1 SST cells with MCs

[p,table,stats] = anovan(Response{:,6},{findgroups(Response(:,2)), ...
    Response{:,15}},'model','interaction','varnames',{'TC','Var15'});
display(Response.Properties.VariableNames{6})                              % amplitude      

[p,table,stats] = anovan(Response{:,7},{findgroups(Response(:,2)), ...
    Response{:,15}},'model','interaction','varnames',{'TC','Var15'});
display(Response.Properties.VariableNames{7})                              % time to peak

%% Making plots for figure

[uniqueRows, ~, groupIndices] = unique([findgroups(Response(:,2)), ...
    Response{:,15}],'rows');

figure
boxchart(groupIndices,Response{:,6})
xticks(unique(groupIndices))
xticklabels({'POM-MC','POM-MMC','VPM-MC','VPM-NMC'})

figure
boxchart(Response{:,15},Response{:,7})

[p,h, stats] = ranksum(Response{Response{:,15}==1,7},Response{Response{:,15}==2,7})


