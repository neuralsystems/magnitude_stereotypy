function pairwiseCosine = computeindividualpaircosinesimilarity(spikeData)
%%COMPUTEINDIVIDUALPAIRCOSINESIMILARITY Calculates the Spearman's rank
%correlation between each individual pair using the vector for all odors.
%
% Usage:
%   pairwiseCosine = COMPUTEINDIVIDUALPAIRCOSINESIMILARITY(spikeData)
%
% Inputs:
%      spikeData: double matrix of size nIndividual x nOdor
%
% Output:
%   pairwiseCosine: correlation of odor vectors between each possible individual pair

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

nInd = size(spikeData, 1);
% find all individual pairs
pairInd = nchoosek(1:nInd, 2);
nPairInd = size(pairInd, 1);
% calculate stereotypy for all individual pairs
pairwiseCosine = zeros(nPairInd, 1);
for iPairInd = 1:nPairInd
    pairwiseCosine(iPairInd) = 1 - pdist([spikeData(pairInd(iPairInd, 1), :); spikeData(pairInd(iPairInd, 2), :)], 'Cosine');
end % for iPairInd
end % computeindividualpaircosinesimilarity