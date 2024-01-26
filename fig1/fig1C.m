data_path = '/Users/celia/Documents/two_context_data';
% CMG
% Most of the legwork is going to be in illustrator. This code will spit
% out the dig percentages in each corner with standard errors.

data = readtable(fullfile(data_path, 'behavior.xlsx'));

animals = unique(data.animal);

%% Collect behavior
behavByDay = cell(length(animals), 3);
c1ByDay = cell(length(animals), 3);
c2ByDay = cell(length(animals),3);
for s = 1:3
    c1Dig = nan(length(animals), 4);
    c2Dig = nan(length(animals), 4);
    allDigs = nan(length(animals),4);
    d = data.session == s;
    if s == 1
        trials2analyze = [9 10 11 12];
    else
        trials2analyze = 1:12;
    end

    for a = 1:length(animals)
        t = data(d & strcmp(data.animal, animals{a}),:);
        context = t.context(trials2analyze);
        dig = categorical(t.firstDig(trials2analyze), {'Corr', 'Geo', 'Feat', 'Wrong'});
        tmpDig = countcats(dig);
        allDigs(a,:) = tmpDig / sum(tmpDig);
        tmpdig = countcats(dig(context == 1))';
        c1Dig(a,:) = tmpdig / sum(tmpdig);
        tmpdig = countcats(dig(context == 2))';
        c2Dig(a,:) = tmpdig / sum(tmpdig);

        behavByDay{a,s} = allDigs(a,:);
        c1ByDay{a,s} = c1Dig(a,:);
        c2ByDay{a,s} = c2Dig(a,:);
    end
end
%% FIG 1C
%
for d = 1:3
    c1Dig = cell2mat(c1ByDay(:,d));
    c2Dig = cell2mat(c2ByDay(:,d));
    c1Sem = std(c1Dig) ./ sqrt(length(c1Dig));
    c2Sem = std(c2Dig) ./ sqrt(length(c2Dig));
    c1Mean = mean(c1Dig);
    c2Mean = mean(c2Dig);

    tmp = [c1Mean; c1Sem];
    fprintf('Day %d:\n\t Context 1 CGNF: %.3f (%.3f), %.3f (%.3f), %.3f (%.3f), %.3f (%.3f)\n', d, tmp(:));
    tmp = [c2Mean; c2Sem];
    fprintf('\t Context 2 CGNF: %.3f (%.3f), %.3f (%.3f), %.3f (%.3f), %.3f (%.3f)\n', tmp(:))
end

