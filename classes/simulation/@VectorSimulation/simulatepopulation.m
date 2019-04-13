function output = simulatepopulation(obj, iInput)
%%SIMULATEPOPULATION Simulates the population of individuals for the
%current stimulus
%
% Usage:
%   SIMULATEPOPULATION(obj, iInput)
%
% Input:
%   iInput: structure specifying the current stimulus
%
% Output:
%   output: cell matrix containing the simulation output for all individuals (size: nIndividual x 1)

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

output = cell(obj.simulationParams.nIndividual, 1);
for iIndividual = 1:obj.simulationParams.nIndividual
    % first resets the input state of all input neurons and then simulates
    % the network for each individual.
    iInput.iIndividual = iIndividual;
    for iInputLayer = obj.population(iIndividual).idInputLayer
        resetinputstate(obj.population(iIndividual).neuronLayer.(iInputLayer{:}), obj.idInputResponse(iIndividual).(iInputLayer{:})(:, iInput.(iInputLayer{:})), obj.currentSeed, iInput);
    end % for iInputNeuron
    output{iIndividual} = simulateneuronnetwork(obj.population(iIndividual));
end % iIndividual
end % simulatepopulation