function pairwiseCorrelation = computeindividualpairpearsoncorrelation(spikeData)
%%COMPUTEINDIVIDUALPAIRPEARSONCORRELATION Calculates the Pearson's
%correlation between each individual pair using the vector for all odors.
%
% Usage:
%   pairwiseCorrelation = COMPUTEINDIVIDUALPAIRPEARSONCORRELATION(spikeData)
%
% Inputs:
%      spikeData: double matrix of size nIndividual x nOdor
%
% Output:
%   pairwiseCorrelation: correlation of odor vectors between each possible individual pair

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

nInd = size(spikeData, 1);
% find all individual pairs
pairInd = nchoosek(1:nInd, 2);
nPairInd = size(pairInd, 1);
% calculate stereotypy for all individual pairs
pairwiseCorrelation = zeros(nPairInd, 1);
for iPairInd = 1:nPairInd
    pairwiseCorrelation(iPairInd) = corr(spikeData(pairInd(iPairInd, 1), :).', spikeData(pairInd(iPairInd, 2), :).', 'Type', 'Pearson', 'rows', 'complete');
end % for iPairInd
end % computeindividualpairpearsoncorrelation