function pairwiseCorrelation = computeindividualpairspearmancorrelation(spikeData)
%%COMPUTEINDIVIDUALPAIRSPEARMANCORRELATION Calculates the Spearman's rank
%correlation between each individual pair using the vector for all odors.
%
% Usage:
%   pairwiseCorrelation = COMPUTEINDIVIDUALPAIRSPEARMANCORRELATION(spikeData)
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
    pairwiseCorrelation(iPairInd) = corr(spikeData(pairInd(iPairInd, 1), :).', spikeData(pairInd(iPairInd, 2), :).', 'Type', 'Spearman');
end % for iPairInd
end % computeindividualpairspearmancorrelation