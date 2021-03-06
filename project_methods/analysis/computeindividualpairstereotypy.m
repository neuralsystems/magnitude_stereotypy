function distanceStereotypy = computeindividualpairstereotypy(spikeData, zeroNaNs)
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
%       zeroNaNs: switch to set nans to 0 or not (Default: true)
%
% Output:
%   distanceStereotypy: distance based stereotypy for each odor and individual pair

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

if nargin < 2
    zeroNaNs = true;
end
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
        distanceStereotypy(iPairInd, iPairOdor) = computepairstereotypy(spikeData(pairInd(iPairInd, :), pairOdor(iPairOdor, :)), zeroNaNs);
    end % for iPairOdor
end % for iPairInd
end % computeindividualpairstereotypy

function stereotypy = computepairstereotypy(spikeData, zeroNaNs)
% calculate squared error between odors for the same individuals
D1 = sum(diff(spikeData) .^ 2);
% calculate squared error between odors for the different individuals
newData = spikeData;
newData(2, :) = newData(2, 2:-1:1);
D2 = sum(diff(newData) .^ 2);
% calculate stereotypy
tempStereotypy = (D2 - D1) ./ (D2 + D1);
% % convert NaNs to 0 (NaNs indicate that both the distances are 0. In such a
% % case we cannot distinguish reliably between odors and individuals and
% % stereotypy should be 0)
if zeroNaNs
    tempStereotypy(isnan(tempStereotypy)) = 0;
end
% report mean stereotypy
stereotypy = tempStereotypy;
end % calculatepairstereotypy