classdef VectorOlfactorySimulation < VectorSimulation
    %%VECTOROLFACTORYSIMULATION Creates an object that extends the VectorSimulation
    %class to simulate the mushroom body for only the olfactory stimulus
    %
    % Usage:
    %   newInstance = VECTOROLFACTORYSIMULATION(params, currentSeed)
    %
    % Inputs:
    %        params: structure containing networkParams and simulationParams that are set during initiation of a simulation instance
    %   currentSeed: scalar specifying the random seed for the current simulation
    %
    % Properties (Required):
    %        currentSeed: scalar specifying the random seed for the current simulation
    %      networkParams: structure specifying all network parameters
    %   simulationParams: structure specifying all simulation parameters
    %              nOdor: scalar specifying the number of odors for which the network can be simulated
    %
    % Properties (Calculated):
    %        population: vector of VectorNeuronNetwork objects (called individuals) that are part of the current simulation
    %   idInputResponse: structure specifying the response matrices for all input neurons to their respective stimulus
    %    simulationData: structure specifying the formatted simulation data for each neuron layer
    %       rateSpiking: structure specifying the spiking rate for each neuron layer
    %
    % Methods:
    %           simulatepopulation: simulates the population of NeuronNetwork objects for all odors
    %           computespikingrate: computes the spiking rate for all the neuron layers
    %    computelayertotalresponse: computes the total response for a neuron layer
    %       computetotallayerinput: computes the input received by a neuron layer
    %   plotsimulationparamsreport: method to plots the report of parameters for the current simulation

    %**********************************************************************%
    % Author: Aarush Mohit Mittal
    % Contact: aarush (dot) mohit (at) gmail (dot) com
    %**********************************************************************%
    
    properties (SetAccess = immutable)
        nOdor
    end % immutable properties
    %**********************************************************************%
    properties (SetAccess = protected)
        simulationData
        rateSpiking
    end % protected properties
    %**********************************************************************%
    properties (Hidden, Access = protected)
        simulationOutput
    end % hidden private properties
    %**********************************************************************%
    methods
        function obj = VectorOlfactorySimulation(params, currentSeed)
            %%VECTOROLFACTORYSIMULATION Class constructor
            obj = obj@VectorSimulation(params, currentSeed);
            obj.nOdor = params.simulationParams.input.(params.networkParams.idInputLayer{1}).nInput;
            simulatepopulation(obj)
            computespikingrate(obj)
        end % class constructor
        %**********************************************************************%
        % declare the implementation of public methods
        simulatepopulation(obj)
        computespikingrate(obj)
        totalResponse = computelayertotalresponse(obj, idLayer)
        totalInputData = computetotallayerinput(obj, idInputLayer, idLayerConnection)
        handleFigure = plotsimulationparamsreport(obj, defaultParams, nSeed)
    end % methods
    %**********************************************************************%
end % VectorOlfactorySimulation