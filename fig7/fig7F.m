%% Side prediction using rates and SVMs
% CMG
% 7/15/2022
data_path = '/Users/celia/Documents/two_context_data'; % change this
load(fullfile(data_path, 'ratemaps-2021_12_1_12_34_49_K1_d4.mat'));

minTrial = 4;
animals = extractfield(ratemaps.data, 'animal');

cPred = cell(length(animals),3);

northDigs = {'Corr'};
southDigs = {'Geo'};

fid = fopen('Fig7Fdata.csv', 'w');
fprintf(fid, 'animal, day, acc\n');

for a = 1:length(animals)
   
    for s = 1:3
        
        if ismember(ratemaps.data(a).animal, 'HGY1_recut')
            if s == 2
                fprintf('HGY1_recut does not have day 2 data\n');
                continue
            elseif s == 3
                fprintf('HGY1_recut does not have day 3 data\n')
                continue
            end
        end
        allcells = findUniqueCells(ratemaps.data(a).session(s));
        contexts = extractfield(ratemaps.data(a).session(s).trial, 'context');
        digs = extractfield(ratemaps.data(a).session(s).trial, 'dig');


        northTrials = find(ismember(digs, northDigs));
        southTrials = find(ismember(digs, southDigs));

        validTrials = [northTrials, southTrials];

        rv = nan(length(validTrials), length(allcells));
        for c = 1:length(allcells)
            mfr = calculateMfrVector(ratemaps.data(a).session(s), allcells{c});
            for tr = 1:length(validTrials)
                wTr = validTrials(tr);
                rv(tr,c) = mfr(wTr);
            end
        end
        allcells = allcells(any(rv,1));
        rv = rv(:, any(rv,1));
        rv(isnan(rv)) = 0; 
        weirdTrial = ~any(rv,2);
        if any(weirdTrial)
            fprintf('%s day %d has a weird trial\n', ratemaps.data(a).animal, s);
        end

        rv = rv(~weirdTrial,:);
        validTrials = validTrials(~weirdTrial);
        

        vDigs = digs(validTrials);
        digClass = cell(1, length(vDigs));
        digClass(ismember(vDigs, northDigs)) = {'North'};
        digClass(ismember(vDigs, southDigs)) = {'South'};
        digClass = categorical(digClass);
        
       mdl = fitcsvm(rv, digClass', 'KernelFunction', 'linear', 'PredictorNames', allcells, 'Prior', 'empirical', 'Leaveout', 'on', 'Verbose',0); 
        
        err = kfoldLoss(mdl, 'mode', 'individual');
        predAcc = 1 - err;
        cPred{a, s} = nanmean(predAcc);
        fprintf(fid, '%s, %d, %.3f\n', animals{a}, s, nanmean(predAcc));
    
    end 
end

avg = nan(3,1);
sem = nan(3,1);
fprintf('Heading Prediction using Support Vector Machines\n')
for d = 1:3
    fprintf('\tDay %d\n', d)
    cdata = cell2mat(cPred(:,d));
    avg(d) = mean(cdata);
    sem(d) = std(cdata) ./ sqrt(length(cdata));
  
    fprintf('\t\t avg = %.4f, sem = %.4f, using %d animals\n',  avg(d), sem(d), length(cdata));
     [p,h] = signrank(cdata, .5, 'tail', 'right');
    fprintf('\t\t\t sign rank test h = %d, p = %.3f\n', h, p);
    [h,p,ci,stats] = ttest(cdata, .5, 'tail', 'right');
    fprintf('\t\t one sample ttest: t(%d) = %.3f, p = %.5f\n', stats.df, stats.tstat, p);
  
end


% Example data 
figure;
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


data1 = cell2mat(cPred(:,1));
data2 = cell2mat(cPred(:,2));
data3 = cell2mat(cPred(:,3));

scatter(repmat(1, length(data1),1), data1, 'ro', 'xjitter', 'rand')
scatter(repmat(2, length(data2),1), data2, 'bo', 'xjitter', 'rand')
scatter(repmat(3, length(data3),1), data3, 'co', 'xjitter', 'rand')

hold off
ylim([0 1])
title('Heading prediction using firing rate');

fclose(fid);