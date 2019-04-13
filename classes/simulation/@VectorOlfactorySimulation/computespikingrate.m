function computespikingrate(obj)
%%COMPUTESPIKINGRATE Computes the spiking rate for all spiking neuron
%layers in the neuron network. Modifies the rateSpiking property of the
%VectorOlfactorySimulation Object
%
% Usage:
%   COMPUTESPIKINGRATE(obj)

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

nIndividual = obj.simulationParams.nIndividual;
obj.rateSpiking = [];
for iLayer = obj.population(1).idNeuronLayer
    obj.rateSpiking.raw.(iLayer{:}) = cell(nIndividual, obj.nOdor);
    obj.rateSpiking.mean.(iLayer{:}).num = 0;
    obj.rateSpiking.mean.(iLayer{:}).all = 0;
    obj.rateSpiking.mean.(iLayer{:}).active = 0;
    obj.rateSpiking.sum.(iLayer{:}).active = 0;
    obj.rateSpiking.mean.(iLayer{:}).baseline = 0;
end
sizeLoop = [nIndividual, obj.nOdor];
nLoop = prod(sizeLoop);
nInput = obj.nOdor * nIndividual;
% calculates the id of active neurons, and spiking rates for active,
% baseline and all neurons
for iLayer = obj.population(1).idNeuronLayer
    for iLoop = 1:nLoop
        [iIndividual, iOdor] = ind2sub(sizeLoop, iLoop);
        rateAll = obj.simulationData.(iLayer{:}){iIndividual, iOdor};
        obj.rateSpiking.raw.(iLayer{:}){iIndividual, iOdor}.allNeuron = rateAll;
        if any(strcmp(iLayer{:}, obj.population(iIndividual).idInputLayer))
            idActive = rateAll > obj.networkParams.neuronLayer.(iLayer{:}).rateBaselineResponse;
            obj.rateSpiking.raw.(iLayer{:}){iIndividual, iOdor}.idActive = idActive;
            obj.rateSpiking.raw.(iLayer{:}){iIndividual, iOdor}.activeNeuron = rateAll(idActive);
            obj.rateSpiking.raw.(iLayer{:}){iIndividual, iOdor}.baselineNeuron = rateAll(~idActive);
        else
            idActive = rateAll > 0;
            obj.rateSpiking.raw.(iLayer{:}){iIndividual, iOdor}.idActive = idActive;
            obj.rateSpiking.raw.(iLayer{:}){iIndividual, iOdor}.activeNeuron = rateAll(idActive);
            obj.rateSpiking.raw.(iLayer{:}){iIndividual, iOdor}.baselineNeuron = rateAll(~idActive);
        end % if iLayer
    end % for iLoop
    obj.rateSpiking.mean.(iLayer{:}).num = arrayfun(@(x) sum(obj.rateSpiking.raw.(iLayer{:}){x}.idActive), 1:nInput);
    obj.rateSpiking.mean.(iLayer{:}).all = arrayfun(@(x) computeultimatemean(obj.rateSpiking.raw.(iLayer{:}){x}.allNeuron), 1:nInput);
    obj.rateSpiking.mean.(iLayer{:}).active = arrayfun(@(x) computeultimatemean(obj.rateSpiking.raw.(iLayer{:}){x}.activeNeuron), 1:nInput);
    obj.rateSpiking.sum.(iLayer{:}).active = arrayfun(@(x) sum(reshape(obj.rateSpiking.raw.(iLayer{:}){x}.activeNeuron, 1, [])), 1:nInput);
    obj.rateSpiking.mean.(iLayer{:}).baseline = arrayfun(@(x) computeultimatemean(obj.rateSpiking.raw.(iLayer{:}){x}.baselineNeuron), 1:nInput);
end % for iLayer
end % computespikingrate