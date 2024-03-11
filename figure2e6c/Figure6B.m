%% Figure 6B

data_path = '/Users/celia/Documents/two_context_data';
calin = load(fullfile(data_path, 'calcium_20220516_170618/analysis_input.mat'));
maps = reformatTbl(calin.analysisInput.MapsData);
calout = load(fullfile(data_path, 'calcium_20220516_170618/analysis_results.mat'));
ba = reformatTbl(calout.analysisResults.BestAligned);
ficells = ba.animalSessionCellName(ba.isStable == 1);
fscells = ba.animalSessionCellName(ba.isStable == 0);

panels = [1,2;1,3; 2, 3]; % 

%% 
clc
for p = 1:size(panels,1)
    fprintf('Day %d to day %d\n', panels(p,1), panels(p,2));
    % find cells registered across days
    cellnames = intersect(maps.registeredCellName(maps.dayUsed == panels(p,1)), ...
        maps.registeredCellName(maps.dayUsed == panels(p,2)));

    cellIdentity = nan(length(cellnames),2);
    for c = 1:length(cellnames)
        cellidx = cellnames{c};
        animalSessionCellName = unique(maps.animalSessionCellName(ismember(maps.registeredCellName, cellidx) & maps.dayUsed == panels(p,1)));
        cellIdentity(c,1) = ismember(animalSessionCellName, ficells);
        animalSessionCellName = unique(maps.animalSessionCellName(ismember(maps.registeredCellName, cellidx) & maps.dayUsed == panels(p,2)));
        cellIdentity(c,2) = ismember(animalSessionCellName, ficells);
    end

    init_fi = cellIdentity(cellIdentity(:,1) == 1,:);
    init_fiNum = size(init_fi,1);
    stayedFI = sum(init_fi(:,2) == 1);
    becameFS = sum(init_fi(:,2) == 0);
    fprintf('\t Out of %d FI cells, %d remained FI and %d became FS\n', init_fiNum, stayedFI, becameFS);
    fprintf('\t\t %.2f%% remained FI, %.2f%% became FS\n', 100*stayedFI/init_fiNum, 100*becameFS/init_fiNum);

    init_fs = cellIdentity(cellIdentity(:,1) == 0,:);
    init_fsNum = size(init_fs,1);
    stayedFS = sum(init_fs(:,2) == 0);
    becameFI = sum(init_fs(:,2) == 1);
    fprintf('\t Out of %d FS cells, %d remained FS and %d became FI\n', init_fsNum, stayedFS, becameFI);
    fprintf('\t\t %.2f%% remained FS, %.2f%% became FI\n', 100*stayedFS/init_fsNum, 100*becameFI/init_fsNum);
end

  