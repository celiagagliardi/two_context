
outputFolder = 'output';
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end
fid = fopen(fullfile('output', 'coherency_table.csv'), 'w');
fprintf(fid, '%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n', 'animal', 'day', 'contextA', 'contextB', 'trialA', 'trialB', 'digA', 'digB', 'animalSessionCellName', 'isStable', 'bmr');
formatSpec = '%s, %d, %d, %d, %d, %d, %s, %s, %s, %d, %d\n';

%% set params and load stuff
data_path = '/Users/celia/Documents/two_context_data'; % change this

calin = load(fullfile(data_path, 'calcium_20220516_170618/analysis_input.mat'));
calout = load(fullfile(data_path, 'calcium_20220516_170618/analysis_results.mat'));
cal_inputT = calin.analysisInput.MapsSquareData;
maps = reformatTbl(cal_inputT);
tmpba = calout.analysisResults.BestAligned;
ba = reformatTbl(tmpba);

minCell = 0; minCorr = -inf;
       
%%
animals = unique(ba.animalName);
outputFolder = 'output';
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end

for a = 1:length(animals)
    atbl = ba(ismember(ba.animalName, animals{a}),:);
    for d = 1:3
        dtbl = atbl(atbl.dayUsed == d,:);
        cellnames = dtbl.animalSessionCellName;
        mtbl = maps(ismember(maps.animalName, animals{a}) & maps.dayUsed == d,:);
        trialId = unique(mtbl.trialId);
        trialPairs = nchoosek(trialId,2);

        contextChanges = [mod(trialPairs(:,1),2), mod(trialPairs(:,2),2)];
        %withinContext = contextChanges(:,1)== contextChanges(:,2);

        %pwRot = nan(size(trialPairs,1),1); % all cells pairwise rotation for each pairwise comparison

        for p = 1:size(trialPairs,1)

            bwRot = nan(length(cellnames),1); % pairwise rotation for current iteration only
            contextA = unique(mtbl.contextId(mtbl.trialId == trialPairs(p,1)));
            contextB = unique(mtbl.contextId(mtbl.trialId == trialPairs(p,2)));
            digA = mtbl.dig{mtbl.trialId == trialPairs(p,1)};
            digB = mtbl.dig{mtbl.trialId == trialPairs(p,2)};

            for c = 1:length(cellnames)
                cellStable = dtbl.isStable(ismember(dtbl.animalSessionCellName, cellnames{c}));


                mapA = mtbl.map{ismember(mtbl.animalSessionCellName, cellnames{c}) & mtbl.trialId == trialPairs(p,1),:};
                mapB = mtbl.map{ismember(mtbl.animalSessionCellName, cellnames{c}) & mtbl.trialId == trialPairs(p,2),:};
                if ~any(mapA, 'all') || ~any(mapB, 'all')
                    bmr = nan;
                end

                rotCorr = nan(1,4);
                rotCorr(1) = corr2(mapA, mapB);
                rotCorr(2) = corr2(mapA, rot90(mapB,1));
                rotCorr(3) = corr2(mapA, rot90(mapB,2));
                rotCorr(4) = corr2(mapA, rot90(mapB,3));
                [rotVal,rot] = max(rotCorr);

                if sum(rotVal == rotCorr) > 1
                    fprintf('marc was right\n');
                end
                if rotVal > minCorr
                    bmr = rot;
                else
                    bmr = nan;
                end

               fprintf(fid, formatSpec, animals{a}, d, contextA, contextB, trialPairs(p,1), trialPairs(p,2), digA, digB, cellnames{c}, cellStable, bmr);

            end

        end
    end
end

fclose(fid)
