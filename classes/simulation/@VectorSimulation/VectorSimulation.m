classdef VectorSimulation < handle
    %%VECTORSIMULATION Creates a simulation object that simulates the
    %network for multiple individuals and inputs. Provides methods for
    %result reporting. Inherits from handle so that new instances can be
    %passed by reference
    %
    % Usage:
    %   newInstance = VECTORSIMULATION(params, currentSeed)
    %
    % Inputs:
    %        params: structure specifying networkParams and simulationParams that are set during initiation of a simulation instance
    %   currentSeed: scalar specifying the random seed for the current simulation
    %
    % Properties (Required):
    %        currentSeed: scalar specifying the random seed for the current simulation
    %      networkParams: structure specifying all network parameters
    %   simulationParams: structure specifying all simulation parameters
    %
    % Properties (Calculated):
    %         population: vector of VectorNeuronNetwork objects (called individuals) that are part of the current simulation
    %    idInputResponse: structure specifying the response matrices for all input neurons to their respective stimulus
    %
    % Methods:
    %   simulatepopulation: simulates the population of all individuals in the population
    %     setinputresponse: sets the response for each input neuron layer in the population
    
    %**********************************************************************%
    % Author: Aarush Mohit Mittal
    % Contact: aarush (dot) mohit (at) gmail (dot) com
    %**********************************************************************%
    properties (SetAccess = immutable)
        currentSeed
        simulationParams
        networkParams
        population@VectorNeuronNetwork
    end % immutable properties
    %**********************************************************************%
    properties (SetAccess = protected)
        idInputResponse
    end % protected properties
    %**********************************************************************%
    methods
        function obj = VectorSimulation(params, currentSeed)
            %%VECTORSIMULATION Class constructor
            obj.currentSeed = currentSeed;
            obj.simulationParams = params.simulationParams;
            obj.networkParams = params.networkParams;
            for iIndividual = 1:params.simulationParams.nIndividual
                obj.population(iIndividual) = VectorNeuronNetwork(params.networkParams, iIndividual, obj.currentSeed);
            end % for iIndividual
            setinputresponse(obj);
        end % class constructor
        %**********************************************************************%
        % declare the implementation of public methods
        setinputresponse(obj)
        output = simulatepopulation(obj, iInput)
    end % methods
    %**********************************************************************%
end % VectorSimulation