%% fig 3C
clear all; clc;
data_path = '/Users/celia/Documents/two_context_data';
load(fullfile(data_path, 'tetrodes_20220316_100902_FINAL/analysis_results.mat'));
tmptbl = analysisResults.BestAligned;
load(fullfile(data_path, 'calcium_20220516_170618/analysis_results.mat'));
tmptbl1 = analysisResults.BestAligned;
tbl = reformatTbl([tmptbl; tmptbl1]);

tbl2 = readtable(fullfile(data_path, 'calcium_tetrode_pairwise_correlation_scatters_by_bestaligned_data_20230614_140623.xlsx'));
tbl2cells = tbl2.Var3;
tbl1cells = tbl.animalSessionCellName;

[c, ia, ib] = intersect(tbl1cells, tbl2cells);

[c, ia] = setdiff(tbl1cells, tbl2cells);
%%

difftbl = tbl(ismember(tbl.animalSessionCellName, c),:);
ctbl = difftbl(end-1,:);
cellname = ctbl.animalSessionCellName;
load(fullfile(data_path, 'calcium_20220516_170618/analysis_input.mat'));
mapsdata = reformatTbl(analysisInput.MapsData);
cmaps = mapsdata(ismember(mapsdata.animalSessionCellName, cellname),:);
figure
for t = 1:12
    tmap = cmaps(cmaps.trialId == t,:);
    if isempty(tmap)
        continue
    end

    subplot(2,6,t)
    pcolor(cell2mat(tmap.map));
    axis ij; daspect([1 1 1]); colormap jet; shading interp;
    title(sprintf('Trial %d', t));
end

%% plot

binsize = 0.02;
figure;
histogram(tbl.bestCorrelation, 'binwidth', binsize, 'normalization', 'probability');
grid on;
xline(0.3);
xlabel('Context Similarity');
ylabel('Probability');

%%

calcells = reformatTbl(tmptbl1);
tetcells = reformatTbl(tmptbl);