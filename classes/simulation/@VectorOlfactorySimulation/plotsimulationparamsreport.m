function handleFigure = plotsimulationparamsreport(obj, defaultParams, nSeed)
%%PLOTSIMULATIONPARAMSREPORT Plots a report of the required simulation
%parameters and returns a handle for the figure
%
% Usage:
%   handleFigure = PLOTSIMULATIONPARAMSREPORT(obj, defaultParams, nSeed)
%
% Inputs:
%   defaultParams: struct containing all the default parameters for the simulation
%           nSeed: number of iterations for each simulation
%
% Output:
%   handleFigure: handle obj to the plotted figure

%**********************************************************************%
% Programmed and Copyright: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

% set initial position of text
global YPos
YPos = 0.5;
% make new figure
handleFigure = figure();
handleAxis = axes(handleFigure);
defaultTextProperties = {'FontSize', 7};
%%----------------------------------------------------------%%
%%%% general parameters
text(handleAxis, 0.1, getY(), sprintf('General Parameters:'), 'FontWeight', 'bold', 'FontSize', 8)
%%% date
text(handleAxis, 0.1, getY(), sprintf('Date: %s', datetime), defaultTextProperties{:})
%%% number of seeds
text(handleAxis, 0.1, getY(), sprintf('Number of Seeds: %d', nSeed), defaultTextProperties{:})
%%% number of individuals
textProperties = [defaultTextProperties, compareProperty(defaultParams.simulationParams.nIndividual, obj.simulationParams.nIndividual)];
text(handleAxis, 0.1, getY(), sprintf('Number of Individuals: %d (%d)', obj.simulationParams.nIndividual, defaultParams.simulationParams.nIndividual), textProperties{:})
%%% number of odors
textProperties = [defaultTextProperties, compareProperty(defaultParams.simulationParams.input.vopn.nInput, obj.nOdor)];
text(handleAxis, 0.1, getY(), sprintf('Number of Odors: %d (%d)', obj.nOdor, defaultParams.simulationParams.input.vopn.nInput), textProperties{:})
getY();
%%----------------------------------------------------------%%
%%%% neuronal layer parameters
text(handleAxis, 0.1, getY(), sprintf('Neuronal Layer Parameters:'), 'FontWeight', 'bold', 'FontSize', 8)
%%% PN layer
textProperties = [defaultTextProperties, compareProperty(defaultParams.networkParams.neuronLayer.vopn.nNeuron, obj.networkParams.neuronLayer.vopn.nNeuron)];
text(handleAxis, 0.1, getY(), sprintf('PN Number: %d (%d)', obj.networkParams.neuronLayer.vopn.nNeuron, defaultParams.networkParams.neuronLayer.vopn.nNeuron), textProperties{:});
%%% KC layer
y = getY();
textProperties = [defaultTextProperties, compareProperty(defaultParams.networkParams.neuronLayer.vokc.nNeuron, obj.networkParams.neuronLayer.vokc.nNeuron)];
text(handleAxis, 0.1, y, sprintf('KC Number: %d (%d)', obj.networkParams.neuronLayer.vokc.nNeuron, defaultParams.networkParams.neuronLayer.vokc.nNeuron), textProperties{:});
textProperties = [defaultTextProperties, compareProperty(defaultParams.networkParams.neuronLayer.vokc.threshold, obj.networkParams.neuronLayer.vokc.threshold)];
text(handleAxis, 5, y, sprintf('KC Response Threshold: %d (%d)', obj.networkParams.neuronLayer.vokc.threshold, defaultParams.networkParams.neuronLayer.vokc.threshold), textProperties{:});
%%% MBON layer
y = getY();
textProperties = [defaultTextProperties, compareProperty(defaultParams.networkParams.neuronLayer.vmbon.nNeuron, obj.networkParams.neuronLayer.vmbon.nNeuron)];
text(handleAxis, 0.1, y, sprintf('MBON Number: %d (%d)', obj.networkParams.neuronLayer.vmbon.nNeuron, defaultParams.networkParams.neuronLayer.vmbon.nNeuron), textProperties{:});
textProperties = [defaultTextProperties, compareProperty(defaultParams.networkParams.neuronLayer.vmbon.threshold, obj.networkParams.neuronLayer.vmbon.threshold)];
text(handleAxis, 5, y, sprintf('MBON Response Threshold: %d (%d)', obj.networkParams.neuronLayer.vmbon.threshold, defaultParams.networkParams.neuronLayer.vmbon.threshold), textProperties{:});
getY();
%%----------------------------------------------------------%%
%%%% PN response parameters
text(handleAxis, 0.1, getY(), sprintf('PN Response Parameters:'), 'FontWeight', 'bold', 'FontSize', 8)
%%% class
textProperties = [defaultTextProperties, compareProperty(func2str(defaultParams.networkParams.neuronLayer.vopn.class), func2str(obj.networkParams.neuronLayer.vopn.class), 'str')];
text(handleAxis, 0.1, getY(), sprintf('Class: %s (%s)', func2str(obj.networkParams.neuronLayer.vopn.class), func2str(defaultParams.networkParams.neuronLayer.vopn.class)), textProperties{:});
%%% response to odor
textProperties = [defaultTextProperties, compareProperty(defaultParams.simulationParams.input.vopn.pInput, obj.simulationParams.input.vopn.pInput)];
text(handleAxis, 0.1, getY(), sprintf('Response Probability to Odor: %.3f (%.3f)', obj.simulationParams.input.vopn.pInput, defaultParams.simulationParams.input.vopn.pInput), textProperties{:});
%%% response matrix variation
textProperties = [defaultTextProperties, compareProperty([defaultParams.simulationParams.input.vopn.vInput.individual, defaultParams.simulationParams.input.vopn.vInput.input], [obj.simulationParams.input.vopn.vInput.individual, obj.simulationParams.input.vopn.vInput.input])];
text(handleAxis, 0.1, getY(), sprintf('Response Variation: Across Individual - %s (%s), Across Odor - %s (%s)', logical2str(obj.simulationParams.input.vopn.vInput.individual), logical2str(defaultParams.simulationParams.input.vopn.vInput.individual), logical2str(obj.simulationParams.input.vopn.vInput.input), logical2str(defaultParams.simulationParams.input.vopn.vInput.input)), textProperties{:});
%%% response matrix function
textProperties = [defaultTextProperties, compareProperty(defaultParams.simulationParams.input.vopn.fnInput, obj.simulationParams.input.vopn.fnInput, 'str')];
text(handleAxis, 0.1, getY(), sprintf('Response Matrix Function: %s (%s)', obj.simulationParams.input.vopn.fnInput, defaultParams.simulationParams.input.vopn.fnInput), textProperties{:});
%%% first odor active spiking rate (only for partitioned case)
if strcmp(func2str(obj.networkParams.neuronLayer.vopn.class), 'VectorOlfactoryProjectionNeuronLayerPartitioned')
    defaultParams.networkParams.neuronLayer.vopn.rateActiveResponseFirst = [10 30];
    textProperties = [defaultTextProperties, compareProperty(defaultParams.networkParams.neuronLayer.vopn.rateActiveResponseFirst, obj.networkParams.neuronLayer.vopn.rateActiveResponseFirst)];
    text(handleAxis, 0.1, getY(), sprintf('Active Neuron Spiking Rate Range for First Odor: %d-%d (%d-%d)', obj.networkParams.neuronLayer.vopn.rateActiveResponseFirst, defaultParams.networkParams.neuronLayer.vopn.rateActiveResponseFirst), textProperties{:});
end
%%% active spiking rate
textProperties = [defaultTextProperties, compareProperty(defaultParams.networkParams.neuronLayer.vopn.rateActiveResponse, obj.networkParams.neuronLayer.vopn.rateActiveResponse)];
text(handleAxis, 0.1, getY(), sprintf('Active Neuron Spiking Rate Range: %d-%d (%d-%d)', obj.networkParams.neuronLayer.vopn.rateActiveResponse, defaultParams.networkParams.neuronLayer.vopn.rateActiveResponse), textProperties{:});
%%% baseline spiking rate
textProperties = [defaultTextProperties, compareProperty(defaultParams.networkParams.neuronLayer.vopn.rateBaselineResponse, obj.networkParams.neuronLayer.vopn.rateBaselineResponse)];
text(handleAxis, 0.1, getY(), sprintf('Baseline Neuron Spiking Rate: %d (%d)', obj.networkParams.neuronLayer.vopn.rateBaselineResponse, defaultParams.networkParams.neuronLayer.vopn.rateBaselineResponse), textProperties{:});
%%% percent noise across individuals
textProperties = [defaultTextProperties, compareProperty(defaultParams.networkParams.neuronLayer.vopn.percentNoise, obj.networkParams.neuronLayer.vopn.percentNoise)];
text(handleAxis, 0.1, getY(), sprintf('Gaussian Noise SD Across Individuals: %d (%d)', obj.networkParams.neuronLayer.vopn.percentNoise, defaultParams.networkParams.neuronLayer.vopn.percentNoise), textProperties{:});
getY();
%%----------------------------------------------------------%%
text(handleAxis, 0.1, getY(), sprintf('Connectivity Parameters:'), 'FontWeight', 'bold', 'FontSize', 8)
%%%% pn to kc connectivity
text(handleAxis, 0.1, getY(), sprintf('PN -> KC:'), 'FontWeight', 'bold', 'FontSize', 7)
%%% connection probability
textProperties = [defaultTextProperties, compareProperty(defaultParams.networkParams.layerConnection.vopn_vokc.pConnection, obj.networkParams.layerConnection.vopn_vokc.pConnection)];
text(handleAxis, 0.1, getY(), sprintf('Connection Probability Between PN and KC: %.2f (%.2f)', mean(obj.networkParams.layerConnection.vopn_vokc.pConnection), mean(defaultParams.networkParams.layerConnection.vopn_vokc.pConnection)), textProperties{:});
%%% randomness
textProperties = [defaultTextProperties, compareProperty(defaultParams.networkParams.layerConnection.vopn_vokc.vIndividual, obj.networkParams.layerConnection.vopn_vokc.vIndividual)];
text(handleAxis, 0.1, getY(), sprintf('Random Connections Between PN and KC Across Individuals? %s (%s)', logical2str(obj.networkParams.layerConnection.vopn_vokc.vIndividual), logical2str(defaultParams.networkParams.layerConnection.vopn_vokc.vIndividual)), textProperties{:});
%%% connectivity function
textProperties = [defaultTextProperties, compareProperty(defaultParams.networkParams.layerConnection.vopn_vokc.fn, obj.networkParams.layerConnection.vopn_vokc.fn, 'str')];
text(handleAxis, 0.1, getY(), sprintf('PN -> KC Connection Matrix Function: %s (%s)', obj.networkParams.layerConnection.vopn_vokc.fn, defaultParams.networkParams.layerConnection.vopn_vokc.fn), textProperties{:});
%%%% kc to mbon connectivity
text(handleAxis, 0.1, getY(), sprintf('KC -> MBON:'), 'FontWeight', 'bold', 'FontSize', 7)
%%% connection probability
textProperties = [defaultTextProperties, compareProperty(defaultParams.networkParams.layerConnection.vokc_vmbon.pConnection, obj.networkParams.layerConnection.vokc_vmbon.pConnection)];
text(handleAxis, 0.1, getY(), sprintf('Connection Probability Between KC and MBON: %.2f (%.2f)', mean(obj.networkParams.layerConnection.vokc_vmbon.pConnection), mean(defaultParams.networkParams.layerConnection.vokc_vmbon.pConnection)), textProperties{:});
%%% randomness
textProperties = [defaultTextProperties, compareProperty(defaultParams.networkParams.layerConnection.vokc_vmbon.vIndividual, obj.networkParams.layerConnection.vokc_vmbon.vIndividual)];
text(handleAxis, 0.1, getY(), sprintf('Random Connections Between KC and MBON Across Individuals? %s (%s)', logical2str(obj.networkParams.layerConnection.vokc_vmbon.vIndividual), logical2str(defaultParams.networkParams.layerConnection.vokc_vmbon.vIndividual)), textProperties{:});
%%% connectivity function
textProperties = [defaultTextProperties, compareProperty(defaultParams.networkParams.layerConnection.vokc_vmbon.fn, obj.networkParams.layerConnection.vokc_vmbon.fn, 'str')];
text(handleAxis, 0.1, getY(), sprintf('KC -> MBON Connection Matrix Function: %s (%s)', obj.networkParams.layerConnection.vokc_vmbon.fn, defaultParams.networkParams.layerConnection.vokc_vmbon.fn), textProperties{:});
%%----------------------------------------------------------%%
% finalise report properties
set(handleAxis, 'YDir', 'reverse')
axis(handleAxis, 'off')
% Add rectangle around text
rectangle(handleAxis, 'Position', [0, 0, 12, getY()])
end
function y = getY()
%%GETY Changes position for each text line
global YPos
y = YPos;
YPos = YPos + 0.5;
end
function propBold = compareProperty(default, new, type)
%%COMPAREPROPERTY Compares the propert values and makes the text bold if
%the default and current properties do not match
if nargin < 3
    type = 'other';
end
switch type
    case 'str'
        if ~strcmp(default, new)
            propBold = {'FontWeight', 'bold'};
        else
            propBold = {};
        end
    otherwise
        if ~isequal(default, new)
            propBold = {'FontWeight', 'bold'};
        else
            propBold = {};
        end
end
end