%% Predict heading and context using only CG trials 
% averaged per animal and two sample t test for prediction rates of FI and
% FS
% Rate as the predictor
% Predict one cell at a time across days.

clear all;
data_path = '/Users/celia/Documents/two_context_data';
load(fullfile(data_path, 'tetrodes_20220316_100902/analysis_input.mat'));
load(fullfile(data_path, 'tetrodes_20220316_100902/analysis_results.mat'));

output_path = '/Users/celia/Documents/two_context_data/figures/SVMpred_celltype';
if ~isfolder(output_path)
    mkdir(output_path);
end

mapsData = reformatTbl(analysisInput.MapsData);
baData = reformatTbl(analysisResults.BestAligned);

northDigs = {'Corr'};
southDigs = {'Geo'};
% add any parameters here
allanimals = unique(baData.animalName);
FIheading = cell(length(allanimals),3);
FSheading = cell(length(allanimals),3);
FIcontext = cell(length(allanimals),3);
FScontext = cell(length(allanimals),3);

fid = fopen('FigS16CDdata.csv', 'w');
fprintf(fid,'animal, day, celltype, heading_prediction, context_prediction\n');
%% compute
for d = 1:3
    for a = 1:length(allanimals)
        celltype = 1; % 1 = FI, 0 = FS;
        allcells = baData.animalSessionCellName(ismember(baData.animalName, allanimals{a}) & baData.dayUsed == d & baData.isStable == celltype);   
        atbl = mapsData(ismember(mapsData.animalName, allanimals{a}) & ismember(mapsData.animalSessionCellName, allcells),:);
        vtbl = atbl(ismember(atbl.dig, [northDigs, southDigs]),:);
        trials = unique(vtbl.trialId);
                
        digClass = cell(length(trials), 1);
        contextClass = nan(length(trials),1);
        for t = 1:length(trials)
            iTrial = trials(t);
            digClass(t) = unique(vtbl.dig(vtbl.trialId == iTrial));
            contextClass(t) = unique(vtbl.contextId(vtbl.trialId == iTrial));
        end
        digClass = categorical(digClass);
        contextClass = categorical(contextClass);

        rv = nan(length(trials), length(allcells));
        for c = 1:length(allcells)
            ctbl = vtbl(ismember(vtbl.animalSessionCellName, allcells{c}),:);
            celltrials = ctbl.trialId;
            for t = 1:length(celltrials)
                iTrial = celltrials(t);
                tIdx = find(trials == iTrial);
                rv(tIdx, c) = ctbl.mfr(ctbl.trialId == iTrial);
            end

        end
        rv(isnan(rv)) = 0;
        
        if ~isempty(digClass)
            hMdl = fitcsvm(rv, digClass, 'kernelfunction', 'linear', 'prior', 'empirical', 'Leaveout', 'on');
            err = kfoldLoss(hMdl);
            acc = 1 - err;
            FIheading{a,d} = acc;            
        end

        if ~isempty(contextClass)
            cMdl = fitcsvm(rv, contextClass, 'kernelfunction', 'linear', 'prior', 'empirical', 'Leaveout', 'on');
            err = kfoldLoss(cMdl);
            acc = 1 - err;
            FIcontext{a,d} = acc;
        end
        
        fprintf(fid, '%s, %d, FI, %.3f, %.3f\n', allanimals{a},d, FIheading{a,d}, FIcontext{a,d} );

        % FS cells
        celltype = 0; % 1 = FI, 0 = FS;
        allcells = baData.animalSessionCellName(ismember(baData.animalName, allanimals{a}) & baData.dayUsed == d & baData.isStable == celltype);   
        atbl = mapsData(ismember(mapsData.animalName, allanimals{a}) & ismember(mapsData.animalSessionCellName, allcells),:);
        vtbl = atbl(ismember(atbl.dig, [northDigs, southDigs]),:);
        trials = unique(vtbl.trialId);
                
        digClass = cell(length(trials), 1);
        contextClass = nan(length(trials),1);
        for t = 1:length(trials)
            iTrial = trials(t);
            digClass(t) = unique(vtbl.dig(vtbl.trialId == iTrial));
            contextClass(t) = unique(vtbl.contextId(vtbl.trialId == iTrial));
        end
        digClass = categorical(digClass');
        contextClass = categorical(contextClass);

        rv = nan(length(trials), length(allcells));
        for c = 1:length(allcells)
            ctbl = vtbl(ismember(vtbl.animalSessionCellName, allcells{c}),:);
            celltrials = ctbl.trialId;
            for t = 1:length(celltrials)
                iTrial = celltrials(t);
                tIdx = find(trials == iTrial);
                rv(tIdx, c) = ctbl.mfr(ctbl.trialId == iTrial);
            end

        end
        rv(isnan(rv)) = 0;

        if ~isempty(digClass)
            hMdl = fitcsvm(rv, digClass, 'kernelfunction', 'linear', 'prior', 'empirical', 'Leaveout', 'on');
            err = kfoldLoss(hMdl);
            acc = 1 - err;
            FSheading{a,d} = acc;
        end

        if ~isempty(contextClass)
            cMdl = fitcsvm(rv, contextClass, 'kernelfunction', 'linear', 'prior', 'empirical', 'Leaveout', 'on');
            err = kfoldLoss(cMdl);
            acc = 1 - err;
            FScontext{a,d} = acc;
        end
        fprintf(fid, '%s, %d, FS, %.3f, %.3f\n', allanimals{a},d, FSheading{a,d}, FScontext{a,d} );
    end
end

fclose(fid);
%% FS vs FI

avg = nan(2,3);
sem = nan(2,3);
figure
for d = 1:3

    fidata = cell2mat(FIheading(:,d));
    fsdata = cell2mat(FSheading(:,d));
    avg(1,d) = nanmean(fidata);
    avg(2,d) = nanmean(fsdata);
    sem(1,d) = std(fidata, 'omitnan') ./ sqrt(sum(~isnan(fidata)));
    sem(2,d) = std(fsdata, 'omitnan') ./ sqrt(sum(~isnan(fsdata)));
   
end
avg = avg';
sem = sem';
subplot(1,2,1)
b = bar(avg, 'grouped');
hold on
% Find the number of groups and the number of bars in each group
[ngroups, nbars] = size(avg);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, avg(:,i), sem(:,i), 'k', 'linestyle', 'none');
end
yline(.5)
hold off
ylim([0 1])
title('Heading Prediction Accuracy');
legend('FI', 'FS')
xlabel('Day');
ylabel('Prediction Accuracy');
    

avg = nan(2,3);
sem = nan(2,3);
for d = 1:3

    fidata = cell2mat(FIcontext(:,d));
    fsdata = cell2mat(FScontext(:,d));
    avg(1,d) = nanmean(fidata);
    avg(2,d) = nanmean(fsdata);
    sem(1,d) = std(fidata, 'omitnan') ./ sqrt(sum(~isnan(fidata)));
    sem(2,d) = std(fsdata, 'omitnan') ./ sqrt(sum(~isnan(fsdata)));  
end
avg = avg';
sem = sem';

subplot(1,2,2)
b = bar(avg, 'grouped');
hold on
% Find the number of groups and the number of bars in each group
[ngroups, nbars] = size(avg);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, avg(:,i), sem(:,i), 'k', 'linestyle', 'none');
end
yline(.5)
hold off
ylim([0 1])
ylabel('Prediction Accuracy');
legend('FI', 'FS')
xlabel('Day');
title('Context Prediction');


%% make boxplots
figure
subplot(1,2,1);
FIheading(cellfun(@isempty, FIheading)) = {nan};
FSheading(cellfun(@isempty, FSheading)) = {nan};
m = cell(1,2);
m{1} = cell2mat(FIheading); m{2} = cell2mat(FSheading);
boxplotGroup(m, 'PrimaryLabels', {'FI', 'FS'}, 'secondaryLabels', {'Day 1', 'Day 2', 'Day 3'})
ylim([0 1.1]);
yline(0.5, '--');
title('Heading Prediction');


subplot(1,2,2);

FIcontext(cellfun(@isempty, FIcontext)) = {nan};
FScontext(cellfun(@isempty, FScontext)) = {nan};
m = cell(1,2);
m{1} = cell2mat(FIcontext); m{2} = cell2mat(FScontext);
boxplotGroup(m, 'PrimaryLabels', {'FI', 'FS'}, 'secondaryLabels', {'Day 1', 'Day 2', 'Day 3'})
ylim([0 1.1]);
yline(0.5, '--')
title('Context Prediction');



