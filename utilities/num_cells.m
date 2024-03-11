%% Num cells 
data_path = '/Users/celia/Documents/two_context_data';

tetout = load(fullfile(data_path, 'tetrodes_20220316_100902_FINAL/analysis_results.mat'));
calout = load(fullfile(data_path, 'calcium_20220516_170618/analysis_results.mat'));
fig2tet = reformatTbl(tetout.analysisResults.NormalBFO90);
fig2cal = reformatTbl(calout.analysisResults.NormalBFO90);
fprintf('Ephys group\n')
for d = 1:3
    dtbl = fig2tet(fig2tet.dayUsed ==d,:);
    num_animals = length(unique(dtbl.animalName));
    numCells = sum(dtbl.all_num_cells);
    fprintf('\tDay %d\n', d);
    fprintf('\t\t %d animals, %d cells\n', num_animals, numCells);
end

fprintf('cal group\n')
for d = 1:3
    dtbl = fig2cal(fig2cal.dayUsed ==d,:);
    num_animals = length(unique(dtbl.animalName));
    numCells = sum(dtbl.all_num_cells);
    fprintf('\tDay %d\n', d);
    fprintf('\t\t %d animals, %d cells\n', num_animals, numCells);
end