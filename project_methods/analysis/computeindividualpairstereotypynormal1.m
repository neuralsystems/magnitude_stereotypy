function distanceStereotypy = computeindividualpairstereotypynormal1(spikeData)
%%COMPUTEINDIVIDUALPAIRSTEREOTYPY Calculates euclidean distance based 
%stereotypy values for each input for each individual and odor pair. It
%finds average D1 (squared error between odors for the same individuals) and
%average D2 (squared error between odors of different individuals). Then
%calculates (D2-D1)/(D2+D1) as the measure of stereotypy.
%
% Usage:
%   distanceStereotypy = COMPUTEINDIVIDUALPAIRSTEREOTYPY(spikeData)
%
% Inputs:
%      spikeData: double matrix of size nIndividual x nOdor
%
% Output:
%   distanceStereotypy: distance based stereotypy for each odor and individual pair

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%
spikeData = normalize(spikeData);
[nInd, nOdor] = size(spikeData);
% find all individual pairs
pairInd = nchoosek(1:nInd, 2);
nPairInd = size(pairInd, 1);
% find all odor pairs
pairOdor = nchoosek(1:nOdor, 2);
nPairOdor = size(pairOdor, 1);
% calculate stereotypy for all individual and odor pairs
distanceStereotypy = zeros(nPairInd, nPairOdor);
for iPairInd = 1:nPairInd
    for iPairOdor = 1:nPairOdor
        distanceStereotypy(iPairInd, iPairOdor) = computepairstereotypy(spikeData(pairInd(iPairInd, :), pairOdor(iPairOdor, :)));
    end % for iPairOdor
end % for iPairInd
end % computeindividualpairstereotypy

function stereotypy = computepairstereotypy(spikeData)
% calculate squared error between odors for the same individuals
D1 = mean(diff(spikeData) .^ 2);
% calculate squared error between odors for the different individuals
newData = spikeData;
newData(2, :) = newData(2, 2:-1:1);
D2 = mean(diff(newData) .^ 2);
% calculate stereotypy
tempStereotypy = (D2 - D1) ./ (D2 + D1);
% convert NaNs to 0 (NaNs indicate that both the distances are 0. In such a
% case we cannot distinguish reliably between odors and individuals and
% stereotypy should be 0)
tempStereotypy(isnan(tempStereotypy)) = 0;
% report mean stereotypy
stereotypy = tempStereotypy;
end % calculatepairstereotypy

function normData = normalize(data)
normData = (data - repmat(min(data, [], 2), 1, size(data, 2))) ./ (repmat(max(data, [], 2), 1, size(data, 2)) - repmat(min(data, [], 2), 1, size(data, 2)));
end