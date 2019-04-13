classdef VectorNeuronNetwork < handle
    %%VECTORNEURONNETWORK Creates a neuron network for a specific
    %individual by adding specified neuron layers. It can also set up
    %connections between the different layers and simulate the network.
    %Inherits from handle so that the class object can be passed by
    %reference (saves time during simulation)
    %
    % Usage:
    %   newInstance = VECTORNEURONNETWORK(params, idIndividual, currentSeed)
    %
    % Inputs:
    %         params: structure of network properties that are set during initiation of a network instance
    %   idIndividual: scalar specifying the id of the individual for which the network is setup
    %    currentSeed: scalar specifying the random seed for the current network conditions
    %
    % Properties (Required):
    %            idIndividual: scalar specifying the id of the individual for which the network is setup
    %             currentSeed: scalar specifying the random seed for the current network conditions
    %    orderLayerSimulation: cell of strings specifying the order in which each neuron layer is simulated
    %         matrixVariation: scalar specifying the percent variation of connections between PNs and KCs across individuals
    %
    % Properties (Calculated):
    %           idNeuronLayer: cell of strings specifying the ids of all neuron layers in the current network
    %            idInputLayer: cell of strings specifying the ids of input neuron layers in the current network
    %             neuronLayer: structure containing the instances for each neuron layer object in the current network
    %   layerConnectionMatrix: structure containing the connection matrices between each neuron layer in the current network
    %
    % Methods:
    %   simulateneuronnetwork: method to simulate all the layers in the neuron network
        
    %**********************************************************************%
    % Author: Aarush Mohit Mittal
    % Contact: aarush (dot) mohit (at) gmail (dot) com
    %**********************************************************************%
    
    properties (SetAccess = immutable)
        idIndividual
        neuronLayer
        idNeuronLayer
        idInputLayer
        matrixVariation
        orderLayerSimulation
        currentSeed
    end % immutable properties
    properties
        layerConnectionMatrix
    end
    methods
        function obj = VectorNeuronNetwork(params, idIndividual, currentSeed)
            %%VECTORNEURONNETWORK Class constructor
            obj.idIndividual = idIndividual;
            obj.orderLayerSimulation = params.orderLayerSimulation;
            obj.idInputLayer = params.idInputLayer;
            obj.matrixVariation = params.matrixVariation;
            obj.currentSeed = currentSeed;
            % Add each neuron layer as a new instance of its corresponding
            % class
            layerParams = params.neuronLayer;
            obj.idNeuronLayer = fieldnames(layerParams).';
            for iNeuronLayer = obj.idNeuronLayer
                obj.neuronLayer.(iNeuronLayer{:}) = layerParams.(iNeuronLayer{:}).class(layerParams.(iNeuronLayer{:}));
            end
            % generate the layer connection matrices for each neuron layer
            obj.layerConnectionMatrix = generatelayerconnectionmatrix(obj, params.layerConnection, setdiff(obj.idNeuronLayer, obj.idInputLayer), obj.matrixVariation);
        end % class constructor
        %**********************************************************************%
        % declare the implementation of public methods
        rateNeuron = simulateneuronnetwork(obj)
    end % methods
    %**********************************************************************%
    methods (Hidden, Access = protected)
        % declare the implementation of hidden protected methods. see
        % individual method help for description
        layerConnection = generatelayerconnectionmatrix(obj, params, idNeuronLayer, matrixVariation)
    end % hidden protected methods
end % VectorNeuronNetwork