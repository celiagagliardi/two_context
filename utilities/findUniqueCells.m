function [cells] = findUniqueCells(sessionRatemap)
% just a helper function to tell me what cells are present in a given
% session. Ratemap structs at the moment only tell you how many cells are
% in a trial, but if a cell doesn't meet threshold for a given trial, it
% won't exist for that trial. Not fair to the cell to pretend it doesn't
% exist.

cells = [];


for t = 1:length(sessionRatemap.trial)
    cells = [cells, sessionRatemap.trial(t).label];
end

cells = unique(cells);