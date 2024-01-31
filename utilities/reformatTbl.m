function [rmTbl] =reformatTbl(tbl, vars2remove)

%tbl = analysisResults.RateMatrices;
animals = unique(tbl.animalName);
rmTbl = [];

for a = 1:length(animals)
    switch animals{a}
        case 'AK42_CA1'
            sessionName = {'d7', 'd8', 'd9'};
        case 'AK74_CA1'
            sessionName = {'d1', 'd2', 'd3'};
        case 'CMG159_recut'
            sessionName = {'s6', 's7', 's8'};
        case 'HGY1_recut'
            sessionName = {'d1'};
        case 'JJ9_CA1'
            sessionName = {'d1', 'd2', 'd3'};
        case 'MG1'
            sessionName = {'d6', 'd7', 'd11'};
        case 'K1_CA1'
            sessionName = {'d1', 'd2', 'd4'};
        case 'CMG129_CA1'
            sessionName = {'s1', 's2', 's3'}; %[1 2 4];
        case 'CMG154_CA1'
            sessionName = {'s1', 's2', 's3'}; %[1 2 3];
        case 'CMG161_CA1'
            sessionName = {'s1', 's2', 's3'}; %[1 2 3];
        case 'CMG162_CA1'
            sessionName = {'s1', 's4', 's5'}; %[1 4 5];
        case 'CMG169_CA1'
            sessionName = {'s1', 's2', 's3'}; %[1 2 3];
    end
    
    for d = 1:length(sessionName)
        x = strcmp(tbl.animalName, animals{a}) & strcmp(tbl.sessionName, sessionName{d});
        dayUsed = repmat(d, sum(x),1);
        tmptbl = [tbl(x,:), array2table(dayUsed)];
        rmTbl = [rmTbl; tmptbl];
        clear tmptbl;
    end
    
end


if nargin > 1
    rmTbl = removevars(rmTbl, vars2remove);
end
