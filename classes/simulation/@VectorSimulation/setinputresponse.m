function setinputresponse(obj)
%%SETINPUTRESPONSE Sets the response of each input neuron to its
%corresponding stimulus. Generates a matrix for each input-input neuron
%pair. Modifies the idInputResponse property of the Simulation object
%
% Usage:
%   SETINPUTRESPONSE(obj)

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

inputParams = obj.simulationParams.input;
nInput = 1;
for iInputLayer = obj.population(1).idInputLayer
    for iIndividual = 1:obj.simulationParams.nIndividual
        % sets the rng for reproducible output and variation across
        % individuals (if vInput.individual is true)
        newSeed = obj.currentSeed + 1e1 * iIndividual * inputParams.(iInputLayer{:}).vInput.individual + 1e2 * nInput;
        rng(newSeed)
        % generates a different matrix for each input if vInput.input is true
        % else duplicates the matrix for all inputs
        if inputParams.(iInputLayer{:}).vInput.input
            sizeMat = [obj.population(iIndividual).neuronLayer.(iInputLayer{:}).nNeuron, inputParams.(iInputLayer{:}).nInput];
            sizeRepmat = 1;
        else
            sizeMat = [obj.population(iIndividual).neuronLayer.(iInputLayer{:}).nNeuron, 1];
            sizeRepmat = [1, inputParams.(iInputLayer{:}).nInput];
        end % if vInput
        obj.idInputResponse(iIndividual).(iInputLayer{:}) = repmat(generaterandommatrix(sizeMat, inputParams.(iInputLayer{:}).pInput, inputParams.(iInputLayer{:}).fnInput), sizeRepmat);
    end % for iIndividual
    nInput = nInput + 1;
end % for iInputLayer
end % setinputresponse