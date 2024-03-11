function [centerOutAngle] = centerOutRadians(map)
% from -180 to 180;

boxSize = size(map);
centerXY = [boxSize(2)/2, boxSize(1)/2];

% calculate center of mass
M = map;
k = prctile(M, 95, 'all');
if isnan(k)
    centerOutAngle = nan;
    return;
end
bw = zeros(size(M));
bw(M>=k) = M(M>=k);

labeledImage = bwlabel(bw);
measurement = regionprops(labeledImage, bw, 'WeightedCentroid');
if length(measurement) > 1 
    [~,tmpA] = max(M, [], 'all');
    [ya, xa] = ind2sub(boxSize, tmpA);
    eucDiff = nan(1,length(measurement));
    for e = 1:length(measurement)
        eucDiff(e) = sqrt((xa - measurement(e).WeightedCentroid(1))^2 + (ya - measurement(e).WeightedCentroid(2))^2);
    end
    [~,Ind] = min(eucDiff);
    com = measurement(Ind).WeightedCentroid;
elseif length(measurement) == 1
    com = measurement.WeightedCentroid;
else 
    centerOutAngle = nan; 
    return;
end

Y = com(2) - centerXY(2); X = com(1) - centerXY(1);
centerOutAngle = atan2(Y,X);