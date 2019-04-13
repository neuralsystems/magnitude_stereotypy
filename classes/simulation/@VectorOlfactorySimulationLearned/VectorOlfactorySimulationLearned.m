classdef VectorOlfactorySimulationLearned < VectorOlfactorySimulation
    %%VECTOROLFACTORYSIMULATIONLEARNED Creates an object that extends the VectorSimulation
    %class to simulate the mushroom body for only the olfactory stimulus
    %
    % Usage:
    %   newInstance = VECTOROLFACTORYSIMULATIONLEARNED(params, currentSeed)
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
    %          pLearning: scalar specifying the probability of KC-MBON connections that are added or deleted for the 2nd odor
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
        pLearning
    end % immutable properties
    methods
        function obj = VectorOlfactorySimulationLearned(params, currentSeed)
            %%VECTOROLFACTORYSIMULATIONLEARNED Class constructor
            obj = obj@VectorOlfactorySimulation(params, currentSeed);
            obj.pLearning = params.simulationParams.pLearning;
            simulatepopulation(obj)
            computespikingrate(obj)
        end % class constructor
        %**********************************************************************%
        % declare the implementation of public methods
        simulatepopulation(obj)
    end % methods
    %**********************************************************************%
end % VectorOlfactorySimulation