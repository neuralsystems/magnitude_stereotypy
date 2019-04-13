function totalResponse = computelayertotalresponse(obj, idLayer)
%%COMPUTETOTALLAYERRESPONSE Computes the total response of a neuron layer
%by summing up the responses for all the neurons in that layer
%
% Usage:
%   totalResponse = COMPUTETOTALLAYERRESPONSE(obj, idLayer)
%
% Inputs:
%   idLayer: string specifying the id of the layer for which the total response is to be calculated
%
% Output:
%   totalResponse: cell containing the total response for the specified layer (size: nIndividual x nOdor)

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

totalResponse = cellfun(@sum, obj.simulationData.(idLayer));
end % computetotallayerresponse