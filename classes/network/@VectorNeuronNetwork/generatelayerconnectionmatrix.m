function layerConnection = generatelayerconnectionmatrix(obj, params, idNeuronLayer, matrixVariation)
%%GENERATELAYERCONNECTIONMATRIX Generates the connection matrices between
%each pair of connected layers in the network
%
% Usage:
%   layerConnection = GENERATELAYERCONNECTIONMATRIX(obj, params, idNeuronLayer, matrixVariation)
%
% Inputs:
%            params: structure containing all the network parameters required for generation of the connection matrices
%     idNeuronLayer: cell of strings containing the id of neuron layers for which the connection matrices need to be generated
%   matrixVariation: scalar specifying the percent variation of connections between PNs and KCs across individuals
%
% Output:
%   layerConnection: structure containing all the connection matrices

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

% a different seed is used for generation of each new matrix for sufficient
% variation across matrices
seedVariation = 1;
for iNeuron = idNeuronLayer
    for iConnection = obj.neuronLayer.(iNeuron{:}).idInputLayer
        % the connection label is generated as an 'inputId_outputId' string
        stringConnection = [iConnection{:}, '_', iNeuron{:}];
        connectionParams = params.(stringConnection);
        if strcmp('vokc', iNeuron{:}) && ~isnan(matrixVariation) && obj.idIndividual > 1
            newSeed = obj.currentSeed + seedVariation * 1e5 + 1e6;
            rng(newSeed)
            matrixSize = [obj.neuronLayer.(iNeuron{:}).nNeuron, obj.neuronLayer.(iConnection{:}).nNeuron];
            originalMatrix = generaterandommatrix(matrixSize, connectionParams.pConnection(obj.idIndividual), connectionParams.fn);
            numCopy = round(numel(originalMatrix) * (1 - matrixVariation));
            newSeed = obj.currentSeed + seedVariation * 1e5 + 1e6 * connectionParams.vIndividual * obj.idIndividual;
            rng(newSeed)
            layerConnection.(stringConnection) = generaterandommatrix(matrixSize, connectionParams.pConnection(obj.idIndividual), connectionParams.fn);
            indexCopy = randperm(numel(originalMatrix), numCopy);
            layerConnection.(stringConnection)(indexCopy) = originalMatrix(indexCopy);
        else
            % the random number generator is set for variation (if required)
            % across individuals and individual connection matrices
            newSeed = obj.currentSeed + seedVariation * 1e5 + 1e6 * connectionParams.vIndividual * obj.idIndividual;
            rng(newSeed)
            matrixSize = [obj.neuronLayer.(iNeuron{:}).nNeuron, obj.neuronLayer.(iConnection{:}).nNeuron];
            layerConnection.(stringConnection) = generaterandommatrix(matrixSize, connectionParams.pConnection(obj.idIndividual), connectionParams.fn);
        end
        seedVariation = seedVariation + 1;
    end % for iConnection
end % for iNeuron
end % generatelayerconnectionmatrix