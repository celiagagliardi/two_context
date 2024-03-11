data_path = '/Users/celia/Documents/two_context_data';
calin = load(fullfile(data_path, 'calcium_20220516_170618/analysis_results.mat'));
tetin = load(fullfile(data_path, 'tetrodes_20220316_100902_FINAL/analysis_results.mat'));
ba = reformatTbl([tetin.analysisResults.BestAligned; calin.analysisResults.BestAligned]);
%%
tbl = ba(ba.isStable == 0,:);
rotlength1 = cellfun(@length, tbl.rotationSequence1)
rotlength2 = cellfun(@length, tbl.rotationSequence2)
numcells = size(tbl,1);
rotidx = find(rotlength1 == 0 | rotlength2 == 0)

fprintf('Of %d FS cells, %d have place fields in only one context: %.3f\n', numcells, length(rotidx), length(rotidx) / numcells);
%%
tbl = ba(ba.isStable == 1,:);
rotlength1 = cellfun(@length, tbl.rotationSequence1);
rotlength2 = cellfun(@length, tbl.rotationSequence2);
numcells = size(tbl,1);
rotidx = find(rotlength1 == 0 | rotlength2 == 0);


fprintf('Of %d FI cells, %d have place fields in only one context: %.3f\n', sum(ba.isStable == 1), length(rotidx), length(rotidx) / sum(ba.isStable == 1));
