ReadTab = ...
    readtable("D:\Cloud\OneDrive - The University of Western Ontario\Rebutal\Reanalysis\TableForAnalysis.xlsx")

firstSixDigits = cellfun(@(x) regexp(x, '\d{6}', 'match', 'once'), ReadTab.Cell_ID, 'UniformOutput', false);

ReadTab.Celltype = categorical(ReadTab.Celltype);
ReadTab.ThalamicNucleus = categorical(ReadTab.ThalamicNucleus);
ReadTab.IC = categorical(ReadTab.IC);
ReadTab.Layer = categorical(ReadTab.Layer);

VPMTab = ReadTab(string(ReadTab.ThalamicNucleus)=="VPM",:);
POMTab = ReadTab(string(ReadTab.ThalamicNucleus)=="POM",:);

VPMMdl = stepwiseglm(VPMTab(:,[2,5,6,7]), 'Binary_Resp ~ 1',...
     'Distribution', 'binomial', 'Link', 'logit','Criterion','Deviance')
POMMdl = stepwiseglm(POMTab(:,[2,5, 6,7]), 'Binary_Resp ~ 1',...
     'Distribution', 'binomial', 'Link', 'logit','Criterion','Deviance')