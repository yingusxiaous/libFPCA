function [V, W] = cleanUp(Ahat, d)
% This function clusters and outputs the correct basis from multiple
% iterations of the FPCA code. Ahat has columns of basis vectors.

% V are the outut vectors.
% W are the unmatched vectors that we're not sure of.

options = statset( 'MaxIter', 1000);
[indices, centres] = kmeans(Ahat', 2*d, 'options', options);
centres = centres';
centres = bsxfun(@rdivide, centres, sum(centres.^2).^0.5);

innerProducts = centres' * centres;
[~, matches] = min(innerProducts);

% Compute the matching vectors which are negatives of each other.
map = [];
leftOver = [];

for j=1:2*d
    if j == matches(matches(j))
        %map{min(j, matches(j))} = max(j, matches(j));
        map(end+1) = min(j, matches(j));
    else
        leftOver(end+1) = j;
    end
end
final = unique(sort(map));

V = centres(:, final);

% An alternative here: use all available data -- flip the signs -- instead
% of just using the first vector.
% V = [];
% for j=1:length(map)
%     if ~isempty(map{j})
%         centre = mean([-Ahat(:, j) Ahat(:, map{j})],2);
%         V(:, end+1) = centre / norm(centre);
%     end
% end

% Cleanup the remaining vectors which we'll use for the next round.
remaining = false(size(indices));

for j=1:length(leftOver)
    remaining = remaining | indices == leftOver(j);
end

W = Ahat(:, remaining);
end