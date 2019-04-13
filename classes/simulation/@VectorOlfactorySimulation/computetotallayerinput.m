function totalInputData = computetotallayerinput(obj, idInputLayer, idLayerConnection)
%%COMPUTETOTALLAYERINPUT Computes the total input received by each neuron
%in a layer from its input layer
%
% Usage:
%   totalInputData = COMPUTETOTALLAYERINPUT(obj, idInputLayer, idLayerConnection)
%
% Inputs:
%        idInputLayer: string specifying the id of the input layer
%   idLayerConnection: string specifying the id of the layer connection in the format input_output layer
%
% Output:
%   totalInputData: cell with each element containing the total number of input spikes (size: nIndividual x nOdor)

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

nIndividual = obj.simulationParams.nIndividual;
totalInputData = cell(nIndividual, obj.nOdor);
sizeLoop = [nIndividual, obj.nOdor];
nLoop = prod(sizeLoop);
for iLoop = 1:nLoop
    [iIndividual, iOdor] = ind2sub(sizeLoop, iLoop);
    totalInputData{iIndividual, iOdor} = obj.population(iIndividual).layerConnectionMatrix.(idLayerConnection) * obj.simulationData.(idInputLayer){iIndividual, iOdor};
end % for iLoop
end % computetotallayerinput