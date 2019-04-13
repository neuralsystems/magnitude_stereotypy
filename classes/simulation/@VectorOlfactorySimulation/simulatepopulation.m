function simulatepopulation(obj)
%%SIMULATEPOPULATION Simulates the population of individuals for all
%inputs. Modifies the simulationData property of the object
%
% Usage:
%   SIMULATEPOPULATION(obj)

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

obj.simulationOutput = cell(obj.nOdor, 1);
for iOdor = 1:obj.nOdor
    iInput.(obj.networkParams.idInputLayer{1}) = iOdor;
    obj.simulationOutput{iOdor} = simulatepopulation@VectorSimulation(obj, iInput);
end % for iOdor
nIndividual = obj.simulationParams.nIndividual;
sizeLoop = [nIndividual, obj.nOdor];
nLoop = prod(sizeLoop);
for idLayer = obj.population(1).idNeuronLayer
    obj.simulationData.(idLayer{:}) = cell(nIndividual, obj.nOdor);
    for iLoop = 1:nLoop
        [iIndividual, iOdor] = ind2sub(sizeLoop, iLoop);
        obj.simulationData.(idLayer{:}){iIndividual, iOdor} ...
            = obj.simulationOutput{iOdor}{iIndividual}.(idLayer{:});
    end % for iLoop
end % for idLayer
end % simulatepopulation