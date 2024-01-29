function [sizeorder, peripheryorder] = randomize_blocks(sizes, peripheries)
% given a set of center sizes and periphery values,
% outputs a randomized ordering of each possible pair
% (every size is matched to every periphery value, once)

nsizes = numel(sizes);
nperipheries = numel(peripheries);

sizeorder = repelem(sizes, 1, nperipheries);
peripheryorder = repmat(peripheries, 1, nsizes);

bothorders = Shuffle([sizeorder; peripheryorder], 1);
sizeorder = bothorders(1, :);
peripheryorder = bothorders(2, :);