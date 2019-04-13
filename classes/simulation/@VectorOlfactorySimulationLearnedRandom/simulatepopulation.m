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

% modify synapses for half the odors randomly
if obj.simulationParams.rateLearning.rate ~= 0
    for iIndividual = 1:nIndividual
        idModifyOdor = randperm(obj.nOdor, round(obj.nOdor / 2));
        if strcmp(obj.simulationParams.rateLearning.type, 'random')
            rateLearning = randsample([-1 1], round(obj.nOdor / 2), true, [0.5 0.5]) * obj.simulationParams.rateLearning.rate;
        else
            rateLearning = ones(1, round(obj.nOdor / 2)) * obj.simulationParams.rateLearning.rate;
        end
        for iOdor  = 1:length(idModifyOdor)
            idActiveKc = obj.simulationOutput{iOdor}{iIndividual}.vokc.' > 0;
            idConnectedKc = obj.population(iIndividual).layerConnectionMatrix.vokc_vmbon;
            if rateLearning(iOdor) > 0;
                idAdd = find(idActiveKc & ~idConnectedKc);
                numAdd = round(numel(idAdd) * rateLearning(iOdor));
                obj.population(iIndividual).layerConnectionMatrix.vokc_vmbon(idAdd(randperm(numel(idAdd), numAdd))) = true;
            else
                idDel = find(idActiveKc & idConnectedKc);
                numDel = round(numel(idDel) * abs(rateLearning(iOdor)));
                obj.population(iIndividual).layerConnectionMatrix.vokc_vmbon(idDel(randperm(numel(idDel), numDel))) = false;
            end
        end
    end
    for iOdor = 1:obj.nOdor
        iInput.(obj.networkParams.idInputLayer{1}) = iOdor;
        obj.simulationOutput{iOdor} = simulatepopulation@VectorSimulation(obj, iInput);
    end % for iOdor
end

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