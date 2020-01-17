function simulateallnetworkvariations()
%%SIMULATEALLNETWORKVARIATIONS Runs all the required simulations for the
%mbon stereotypy project
%
% Usage:
%   SIMULATEALLNETWORKVARIATIONS()
%
% Simulations Defined:
%   01. Default Network with No Changes
%   02. Default Network with 100 Odors
%   03. Default Network with learning
%   04. Partitioned Network with No Changes
%   05. Partitioned Network Cases with First Odor at [10 30] and Second Odor having Different Firing Range
%   06. Partitioned Network Cases with Both Odors having Same Firing Ranges
%   07. Linear Partitioned Network with First Odor at [10 30] and Second Odor having Different Firing Range
%   08. Partitioned Network Cases with First Odor having 25 Active PNs and Second Odor having Different Number of Active PNs
%   09. Partitioned Network Cases with Both Odors having Same Number of Active PNs
%   10. Linear Partitioned Network Cases with First Odor having 25 Active PNs and Second Odor having Different Number of Active PNs
%   11. Shuffled Network with No Changes
%   12. Default Network with Variation in PN Mean Firing Rate
%   13. Default Network with Variation in KC Number
%   14. Default Network Cases with Variation in KC-MBON Connection Probability
%   15. Default Network Cases with Variation in percent variation in PN-KC connections across individuals and in KC-MBON Connection Probability
%   16. Default Network with 100 odors and Same Connections Across Individuals Without Noise
%   17. Default Network with 100 odors and Same Connections Across Individuals With Noise
%   18. Default Network with 100 odors and Same Connections Across Individuals Without Noise with Odors varying Across Individuals
%   19. Default Network with 100 odors and noise
%   20. Default network with 2 individuals
%   21. Default network with 3 individuals
%   22. Default network with 4 individuals
%   23. Default network with new seed
%   24. Partitioned Network Cases with Both Odors having Different Input Drives
%   25. Default Network with Variation in PN-KC Connection Probability
%   26. Default Network with Variation in KC Response Threshold
%   27. Default Network Cases with Variation in number of PNs
%   28. Default Network Cases with Variation in PN response probability

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

% define common parameters
defaultParams = generatedefaultparams(@defaultparamsvectornetwork, @defaultparamsvectorsimulation);
pathCommon = 'mbon_stereotypy_final\';
%---------------------------------------------------------------------------------------------------%

%%%%%%%%%%--Default Network with No Changes--%%%%%%%%%%
nSeed = 100; % number of network iterations
rng('default'); % reset random number generator for reproducible results
idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
% make simulation folder
pathSimulation = [pathCommon, 'default_network\'];
pathResult = createresultfolder(pathSimulation);
% set new params
newParams = defaultParams;
% simulate network
S = cell(nSeed, 1);
for iSeed = 1:nSeed
    fprintf('Iteration %d: ', iSeed);
    timeStamp = tic;
    S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
    toc(timeStamp)
end % for iSeed
save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
%---------------------------------------------------------------------------------------------------%

%%%%%%%%%%--Default Network with 100 Odors--%%%%%%%%%%
nSeed = 100;
rng('default'); % reset random number generator for reproducible results
idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
% make simulation folder
pathSimulation = [pathCommon, 'default_network_odor_100\'];
pathResult = createresultfolder(pathSimulation);
% set new params
newParams = defaultParams;
newParams.simulationParams.input.vopn.nInput = 100;
% simulate network
for iSeed = 1:nSeed
    fprintf('Iteration %d: ', iSeed);
    timeStamp = tic;
    S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
    toc(timeStamp)
end % for iSeed
save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
%---------------------------------------------------------------------------------------------------%

%%%%%%%%%%--Default Network with learning--%%%%%%%%%%
nSeed = 100;
rng('default'); % reset random number generator for reproducible results
% set new params
newParams = defaultParams;
newParams.simulationParams.input.vopn.nInput = 100;
% set variation ranges
learningRange = 0:0.2:1;
for rateLearning = learningRange
    idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
    % set kc to mbon connection probability of individual 2
    newParams.simulationParams.rateLearning.type = 'random';
    newParams.simulationParams.rateLearning.rate = rateLearning;
    % make simulation folder
    pathSimulation = [pathCommon, sprintf('default_network_odor_100_learn_random_both_p_%.2f\\', rateLearning)];
    pathResult = createresultfolder(pathSimulation);
    % simulate network
    S = cell(nSeed, 1);
    for iSeed = 1:nSeed
        fprintf('Iteration %d: ', iSeed);
        timeStamp = tic;
        S{iSeed} = VectorOlfactorySimulationLearnedRandom(newParams, idSeed(iSeed));
        toc(timeStamp)
    end % for iSeed
    save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
end
%---------------------------------------------------------------------------------------------------%

%%%%%%%%%%--Partitioned Network with No Changes--%%%%%%%%%%
nSeed = 100;
rng('default'); % reset random number generator for reproducible results
idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
% make simulation folder
pathSimulation = [pathCommon, 'partitioned_network\'];
pathResult = createresultfolder(pathSimulation);
% set new params
newParams = defaultParams;
% change PN class to partitioned
newParams.networkParams.neuronLayer.vopn.class = @VectorOlfactoryProjectionNeuronLayerPartitioned;
% change PN firing probability so that each odor has the same number of active PNs (required for partitioning)
newParams.simulationParams.input.vopn.fnInput = 'pseudo-random-column';
% add separate firing rate for the first odor
newParams.networkParams.neuronLayer.vopn.rateActiveResponseFirst = [10 30];
% simulate network
S = cell(nSeed, 1);
for iSeed = 1:nSeed
    fprintf('Iteration %d: ', iSeed);
    timeStamp = tic;
    S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
    toc(timeStamp)
end % for iSeed
save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
%---------------------------------------------------------------------------------------------------%

%%%%%%%%%%--Partitioned Network Cases with First Odor at [10 30] and Second Odor having Different Firing Range--%%%%%%%%%%
nSeed = 100;
rng('default'); % reset random number generator for reproducible results
idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
% set new params
newParams = defaultParams;
% change PN class to partitioned
newParams.networkParams.neuronLayer.vopn.class = @VectorOlfactoryProjectionNeuronLayerPartitioned;
% change PN firing probability so that each odor has the same number of active PNs (required for partitioning)
newParams.simulationParams.input.vopn.fnInput = 'pseudo-random-column';
% add separate firing rate for the first odor
newParams.networkParams.neuronLayer.vopn.rateActiveResponseFirst = [10 30];
% set variation ranges
rateMin = 10;
rateMax = 30:5:80;
for iRate = 1:length(rateMax)
    % set firing range for second odor
    newParams.networkParams.neuronLayer.vopn.rateActiveResponse = [rateMin, rateMax(iRate)];
    % make simulation folder
    pathSimulation = [pathCommon, sprintf('partitioned_network_both_odor_different_firing_range_%d_%d\\', rateMin, rateMax(iRate))];
    pathResult = createresultfolder(pathSimulation);
    % simulate network
    S = cell(nSeed, 1);
    for iSeed = 1:nSeed
        fprintf('Iteration %d: ', iSeed);
        timeStamp = tic;
        S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
        toc(timeStamp)
    end % for iSeed
    save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
end
%---------------------------------------------------------------------------------------------------%

%%%%%%%%%%--Partitioned Network Cases with Both Odors having Same Firing Ranges--%%%%%%%%%%
nSeed = 100;
rng('default'); % reset random number generator for reproducible results
idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
% set new params
newParams = defaultParams;
% change PN class to partitioned
newParams.networkParams.neuronLayer.vopn.class = @VectorOlfactoryProjectionNeuronLayerPartitioned;
% change PN firing probability so that each odor has the same number of active PNs (required for partitioning)
newParams.simulationParams.input.vopn.fnInput = 'pseudo-random-column';
rateMin = 10;
rateMax = 30:5:80;
for iRate = 1:length(rateMax)
    % set firing range for both odors
    newParams.networkParams.neuronLayer.vopn.rateActiveResponseFirst = [rateMin, rateMax(iRate)];
    newParams.networkParams.neuronLayer.vopn.rateActiveResponse = [rateMin, rateMax(iRate)];
    % make simulation folder
    pathSimulation = [pathCommon, sprintf('partitioned_network_both_odor_same_firing_range_%d_%d\\', rateMin, rateMax(iRate))];
    pathResult = createresultfolder(pathSimulation);
    % simulate network
    S = cell(nSeed, 1);
    for iSeed = 1:nSeed
        fprintf('Iteration %d: ', iSeed);
        timeStamp = tic;
        S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
        toc(timeStamp)
    end % for iSeed
    save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
end
%---------------------------------------------------------------------------------------------------%

%%%%%%%%%%--Linear Partitioned Network Cases with First Odor at [10 30] and Second Odor having Different Firing Range--%%%%%%%%%%
nSeed = 100;
rng('default'); % reset random number generator for reproducible results
idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
% set new params
newParams = defaultParams;
% change KCs and MBON class to linear
newParams.networkParams.neuronLayer.vokc.class = @VectorSpikingNeuronLayerLinear;
newParams.networkParams.neuronLayer.vmbon.class = @VectorSpikingNeuronLayerLinear;
% add slope value to the params
newParams.networkParams.neuronLayer.vokc.slope = 1.755;
newParams.networkParams.neuronLayer.vmbon.slope = 1.755;
% change PN class to partitioned
newParams.networkParams.neuronLayer.vopn.class = @VectorOlfactoryProjectionNeuronLayerPartitioned;
% change PN firing probability so that each odor has exact number of active PNs (required for partitioning)
newParams.simulationParams.input.vopn.fnInput = 'pseudo-random-column';
% add separate firing rate for the first odor
newParams.networkParams.neuronLayer.vopn.rateActiveResponseFirst = [10 30];
% set variation ranges
rateMin = 10;
rateMax = 30:5:80;
for iRate = 1:length(rateMax)
    % set firing range for second odor
    newParams.networkParams.neuronLayer.vopn.rateActiveResponse = [rateMin, rateMax(iRate)];
    % make simulation folder
    pathSimulation = [pathCommon, sprintf('linear_partitioned_network_both_odor_different_firing_range_%d_%d\\', rateMin, rateMax(iRate))];
    pathResult = createresultfolder(pathSimulation);
    % simulate network
    S = cell(nSeed, 1);
    for iSeed = 1:nSeed
        fprintf('Iteration %d: ', iSeed);
        timeStamp = tic;
        S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
        toc(timeStamp)
    end % for iSeed
    save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
end
%---------------------------------------------------------------------------------------------------%

%%%%%%%%%%--Partitioned Network Cases with First Odor having 25 Active PNs and Second Odor having Different Number of Active PNs--%%%%%%%%%%
nSeed = 100;
rng('default'); % reset random number generator for reproducible results
idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
% set new params
newParams = defaultParams;
% change PN class to partitioned
newParams.networkParams.neuronLayer.vopn.class = @VectorOlfactoryProjectionNeuronLayerPartitionedNumber;
% change PN firing probability so that each odor has the same number of active PNs (required for partitioning)
newParams.simulationParams.input.vopn.fnInput = 'pseudo-random-column';
nActivePn = 25:45;
for iSim = 1:length(nActivePn)
    % set active pn number for both odors
    newParams.networkParams.neuronLayer.vopn.nActiveFirst = nActivePn(iSim);
    % make simulation folder
    pathSimulation = [pathCommon, sprintf('partitioned_network_both_odor_different_number_%d\\', nActivePn(iSim))];
    pathResult = createresultfolder(pathSimulation);
    % simulate network
    S = cell(nSeed, 1);
    for iSeed = 1:nSeed
        fprintf('Iteration %d: ', iSeed);
        timeStamp = tic;
        S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
        toc(timeStamp)
    end % for iSeed
    save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
end
%---------------------------------------------------------------------------------------------------%

%%%%%%%%%%--Partitioned Network Cases with Both Odors having Same Number of Active PNs--%%%%%%%%%%
nSeed = 100;
rng('default'); % reset random number generator for reproducible results
idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
% set new params
newParams = defaultParams;
% change PN class to partitioned
newParams.networkParams.neuronLayer.vopn.class = @VectorOlfactoryProjectionNeuronLayerPartitioned;
% change PN firing probability so that each odor has the same number of active PNs (required for partitioning)
newParams.simulationParams.input.vopn.fnInput = 'pseudo-random-column';
% set firing rate range for both odors
newParams.networkParams.neuronLayer.vopn.rateActiveResponseFirst = [10 30];
nActivePn = 25:45;
for iSim = 1:length(nActivePn)
    % set active pn probability for both odors
    newParams.simulationParams.input.vopn.pInput = nActivePn(iSim) / 50;
    % make simulation folder
    pathSimulation = [pathCommon, sprintf('partitioned_network_both_odor_same_number_%d\\', nActivePn(iSim))];
    pathResult = createresultfolder(pathSimulation);
    % simulate network
    S = cell(nSeed, 1);
    for iSeed = 1:nSeed
        fprintf('Iteration %d: ', iSeed);
        timeStamp = tic;
        S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
        toc(timeStamp)
    end % for iSeed
    save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
end
%---------------------------------------------------------------------------------------------------%

%%%%%%%%%%--Linear Partitioned Network Cases with First Odor having 25 Active PNs and Second Odor having Different Number of Active PNs--%%%%%%%%%%
nSeed = 100;
rng('default'); % reset random number generator for reproducible results
idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
% set new params
newParams = defaultParams;
% change KCs and MBON class to linear
newParams.networkParams.neuronLayer.vokc.class = @VectorSpikingNeuronLayerLinear;
newParams.networkParams.neuronLayer.vmbon.class = @VectorSpikingNeuronLayerLinear;
% add slope value to the params
newParams.networkParams.neuronLayer.vokc.slope = 1.755;
newParams.networkParams.neuronLayer.vmbon.slope = 1.755;
% change PN class to partitioned
newParams.networkParams.neuronLayer.vopn.class = @VectorOlfactoryProjectionNeuronLayerPartitionedNumber;
% change PN firing probability so that each odor has exact number of active PNs (required for partitioning)
newParams.simulationParams.input.vopn.fnInput = 'pseudo-random-column';
nActivePn = 25:45;
for iSim = 1:length(nActivePn)
    % set active pn number for both odors
    newParams.networkParams.neuronLayer.vopn.nActiveFirst = nActivePn(iSim);
    % make simulation folder
    pathSimulation = [pathCommon, sprintf('linear_partitioned_network_both_odor_different_number_%d\\', nActivePn(iSim))];
    pathResult = createresultfolder(pathSimulation);
    % simulate network
    S = cell(nSeed, 1);
    for iSeed = 1:nSeed
        fprintf('Iteration %d: ', iSeed);
        timeStamp = tic;
        S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
        toc(timeStamp)
    end % for iSeed
    save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
end
%---------------------------------------------------------------------------------------------------%

%%%%%%%%%%--Shuffled Network with No Changes--%%%%%%%%%%
nSeed = 100;
rng(1); % reset random number generator for reproducible results
idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
% make simulation folder
pathSimulation = [pathCommon, 'shuffled_network\'];
pathResult = createresultfolder(pathSimulation);
% set new params
newParams = defaultParams;
% change PN class to shuffled
newParams.networkParams.neuronLayer.vopn.class = @VectorOlfactoryProjectionNeuronLayerShuffled;
% make active PNs same across individuals (made different later)
newParams.simulationParams.input.vopn.vInput.input = false;
% simulate network
S = cell(nSeed, 1);
for iSeed = 1:nSeed
    fprintf('Iteration %d: ', iSeed);
    timeStamp = tic;
    S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
    toc(timeStamp)
end % for iSeed
save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
%---------------------------------------------------------------------------------------------------%

%%%%%%%%%%--Default Network Cases with Variation in PN Mean Firing Rate--%%%%%%%%%%
nSeed = 100;
rng('default'); % reset random number generator for reproducible results
idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
% set new params
newParams = defaultParams;
% set variation ranges
pnMeanVarRange = -10:5:50;
for nVar = pnMeanVarRange
    % set PN firing range
    newParams.networkParams.neuronLayer.vopn.rateActiveResponse = [10 30] + nVar;
    % make simulation folder
    pathSimulation = [pathCommon, sprintf('default_network_pn_firing_rate_mean_%d\\', 20 + nVar)];
    pathResult = createresultfolder(pathSimulation);
    % simulate network
    S = cell(nSeed, 1);
    for iSeed = 1:nSeed
        fprintf('Iteration %d: ', iSeed);
        timeStamp = tic;
        S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
        toc(timeStamp)
    end % for iSeed
    save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
end
%---------------------------------------------------------------------------------------------------%

%%%%%%%%%%--Default Network Cases with Variation in KC Number--%%%%%%%%%%
nSeed = 100;
rng('default'); % reset random number generator for reproducible results
idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
% set new params
newParams = defaultParams;
% set variation ranges
nKcRange = [100 250 500 1e3:1e3:9e3 1e4:2e3:1.8e4];
for nKc = nKcRange
    % set KC number
    newParams.networkParams.neuronLayer.vokc.nNeuron = nKc;
    % make simulation folder
    pathSimulation = [pathCommon, sprintf('default_network_num_kc_%d\\', nKc)];
    pathResult = createresultfolder(pathSimulation);
    % simulate network
    S = cell(nSeed, 1);
    for iSeed = 1:nSeed
        fprintf('Iteration %d: ', iSeed);
        timeStamp = tic;
        S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
        toc(timeStamp)
    end % for iSeed
    save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
end
%---------------------------------------------------------------------------------------------------%

%%%%%%%%%%--Default Network Cases with Variation in KC-MBON Connection Probability--%%%%%%%%%%
nSeed = 100;
rng('default'); % reset random number generator for reproducible results
idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
% set new params
newParams = defaultParams;
% set variation ranges
cKcMbonRange = [0.01, 0.02:0.02:0.08, 0.1:0.1:1];
for c = cKcMbonRange
    % set KC number
    newParams.networkParams.layerConnection.vokc_vmbon.pConnection = [c, c];
    % make simulation folder
    pathSimulation = [pathCommon, sprintf('default_network_c_kc_mbon_%.2f\\', c)];
    pathResult = createresultfolder(pathSimulation);
    % simulate network
    S = cell(nSeed, 1);
    for iSeed = 1:nSeed
        fprintf('Iteration %d: ', iSeed);
        timeStamp = tic;
        S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
        toc(timeStamp)
    end % for iSeed
    save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
end
%---------------------------------------------------------------------------------------------------%

%%%%%%%%%%--Default Network Cases with Variation in percent variation in PN-KC connections across individuals and in KC -> MBON Connection Probability--%%%%%%%%%%
nSeed = 100;
rng('default'); % reset random number generator for reproducible results
idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
% set new params
newParams = defaultParams;
% set variation ranges
variationRange = 10 .^ linspace(-2, 0, 21);
cKcMbonRange = 10 .^ linspace(-2, 0, 21);
for v = variationRange
    for c = cKcMbonRange
        % set KC number
        newParams.networkParams.layerConnection.vokc_vmbon.pConnection = [c, c];
        newParams.networkParams.matrixVariation = v;
        % make simulation folder
        pathSimulation = [pathCommon, sprintf('default_network_c_kc_mbon_%.3f_var_pn_kc_%.3f\\', c, v)];
        pathResult = createresultfolder(pathSimulation);
        % simulate network
        S = cell(nSeed, 1);
        for iSeed = 1:nSeed
            fprintf('Iteration %d: ', iSeed);
            timeStamp = tic;
            S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
            toc(timeStamp)
        end % for iSeed
        idSeed = sort(randperm(1e5, nSeed));
        save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
    end
end
%---------------------------------------------------------------------------------------------------%
%%%%%%%%%%--Default Network with 100 odors and Same Connections Across Individuals Without Noise--%%%%%%%%%%
nSeed = 100;
rng('default'); % reset random number generator for reproducible results
idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
% make simulation folder
pathSimulation = [pathCommon, 'default_network_odor_100_same_connection\'];
pathResult = createresultfolder(pathSimulation);
% set new params
newParams = defaultParams;
newParams.simulationParams.input.vopn.nInput = 100;
% make PN -> KC connections across individuals same
newParams.networkParams.layerConnection.vopn_vokc.vIndividual = false;
% simulate network
S = cell(nSeed, 1);
for iSeed = 1:nSeed
    fprintf('Iteration %d: ', iSeed);
    timeStamp = tic;
    S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
    toc(timeStamp)
end % for iSeed
save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
%---------------------------------------------------------------------------------------------------%

%%%%%%%%%%--Default Network with 100 odors and Same Connections Across Individuals With Noise--%%%%%%%%%%
nSeed = 100;
rng('default'); % reset random number generator for reproducible results
idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
% make simulation folder
pathSimulation = [pathCommon, 'default_network_odor_100_same_connection_with_noise\'];
pathResult = createresultfolder(pathSimulation);
% set new params
newParams = defaultParams;
newParams.simulationParams.input.vopn.nInput = 100;
% make PN -> KC connections across individuals same
newParams.networkParams.layerConnection.vopn_vokc.vIndividual = false;
% add some gaussian noise in the PN responses
newParams.networkParams.neuronLayer.vopn.percentNoise = 2;
% add some gaussian noise in the KC responses
newParams.networkParams.neuronLayer.vokc.class = @VectorSpikingNeuronLayerNoisy;
newParams.networkParams.neuronLayer.vokc.percentNoise = 2;
% add some gaussian noise in the MBON responses
newParams.networkParams.neuronLayer.vmbon.class = @VectorSpikingNeuronLayerNoisy;
newParams.networkParams.neuronLayer.vmbon.percentNoise = 2;
% simulate network
S = cell(nSeed, 1);
for iSeed = 1:nSeed
    fprintf('Iteration %d: ', iSeed);
    timeStamp = tic;
    S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
    toc(timeStamp)
end % for iSeed
save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
%---------------------------------------------------------------------------------------------------%

%%%%%%%%%%--Default Network with 100 odors and Same Connections Across Individuals Without Noise with Odors varying Across Individuals--%%%%%%%%%%
nSeed = 100;
rng('default'); % reset random number generator for reproducible results
idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
% make simulation folder
pathSimulation = [pathCommon, 'default_network_odor_100_same_connection_vary_odor_across_individual\'];
pathResult = createresultfolder(pathSimulation);
% set new params
newParams = defaultParams;
newParams.simulationParams.input.vopn.nInput = 100;
% make PN -> KC connections across individuals same
newParams.networkParams.layerConnection.vopn_vokc.vIndividual = false;
% make odors vary across individuals
newParams.simulationParams.input.vopn.vInput.individual = true;
% simulate network
S = cell(nSeed, 1);
for iSeed = 1:nSeed
    fprintf('Iteration %d: ', iSeed);
    timeStamp = tic;
    S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
    toc(timeStamp)
end % for iSeed
save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
%---------------------------------------------------------------------------------------------------%

%%%%%%%%%%--Default Network with 100 odors and Noise--%%%%%%%%%%
nSeed = 100;
rng('default'); % reset random number generator for reproducible results
% set new params
newParams = defaultParams;
newParams.simulationParams.input.vopn.nInput = 100;
% set noise range
noiseRange = [2 4:4:12];
for noise = noiseRange
    idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
    % make simulation folder
    pathSimulation = [pathCommon, sprintf('default_network_odor_100_with_noise_%d\\', noise)];
    pathResult = createresultfolder(pathSimulation);
    % add some gaussian noise in the PN responses
    newParams.networkParams.neuronLayer.vopn.percentNoise = noise;
    % add some gaussian noise in the KC responses
    newParams.networkParams.neuronLayer.vokc.class = @VectorSpikingNeuronLayerNoisy;
    newParams.networkParams.neuronLayer.vokc.percentNoise = noise;
    % add some gaussian noise in the MBON responses
    newParams.networkParams.neuronLayer.vmbon.class = @VectorSpikingNeuronLayerNoisy;
    newParams.networkParams.neuronLayer.vmbon.percentNoise = noise;
    % simulate network
    S = cell(nSeed, 1);
    for iSeed = 1:nSeed
        fprintf('Iteration %d: ', iSeed);
        timeStamp = tic;
        S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
        toc(timeStamp)
    end % for iSeed
    save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
end
%---------------------------------------------------------------------------------------------------%
%%%%%%%%%%--Default Network with 2 individuals--%%%%%%%%%%
nSeed = 100; % number of network iterations
rng('default'); % reset random number generator for reproducible results
idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
% make simulation folder
pathSimulation = [pathCommon, 'default_network_ind_2\'];
pathResult = createresultfolder(pathSimulation);
% set new params
newParams = defaultParams;
newParams.simulationParams.nIndividual = 2;
newParams.networkParams.layerConnection.vokc_vmbon.pConnection = 0.5 * ones(1, newParams.simulationParams.nIndividual);
newParams.networkParams.layerConnection.vopn_vokc.pConnection = 0.14 * ones(1, newParams.simulationParams.nIndividual);
% simulate network
S = cell(nSeed, 1);
for iSeed = 1:nSeed
    fprintf('Iteration %d: ', iSeed);
    timeStamp = tic;
    S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
    toc(timeStamp)
end % for iSeed
save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
%---------------------------------------------------------------------------------------------------%

%%%%%%%%%%--Default Network with 3 individuals--%%%%%%%%%%
nSeed = 100;
rng('default'); % reset random number generator for reproducible results
idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
% make simulation folder
pathSimulation = [pathCommon, 'default_network_ind_3\'];
pathResult = createresultfolder(pathSimulation);
% set new params
newParams = defaultParams;
newParams.simulationParams.nIndividual = 3;
newParams.networkParams.layerConnection.vokc_vmbon.pConnection = 0.5 * ones(1, newParams.simulationParams.nIndividual);
newParams.networkParams.layerConnection.vopn_vokc.pConnection = 0.14 * ones(1, newParams.simulationParams.nIndividual);
% simulate network
S = cell(nSeed, 1);
for iSeed = 1:nSeed
    fprintf('Iteration %d: ', iSeed);
    timeStamp = tic;
    S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
    toc(timeStamp)
end % for iSeed
save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
%---------------------------------------------------------------------------------------------------%

%%%%%%%%%%--Default Network with 4 individuals--%%%%%%%%%%
nSeed = 100;
rng('default'); % reset random number generator for reproducible results
idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
% make simulation folder
pathSimulation = [pathCommon, 'default_network_ind_4\'];
pathResult = createresultfolder(pathSimulation);
% set new params
newParams = defaultParams;
newParams.simulationParams.nIndividual = 4;
newParams.networkParams.layerConnection.vokc_vmbon.pConnection = 0.5 * ones(1, newParams.simulationParams.nIndividual);
newParams.networkParams.layerConnection.vopn_vokc.pConnection = 0.14 * ones(1, newParams.simulationParams.nIndividual);
% simulate network
S = cell(nSeed, 1);
for iSeed = 1:nSeed
    fprintf('Iteration %d: ', iSeed);
    timeStamp = tic;
    S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
    toc(timeStamp)
end % for iSeed
save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
%---------------------------------------------------------------------------------------------------%

%%%%%%%%%%--Default Network with New Seed--%%%%%%%%%%
nSeed = 100;
rng(100); % reset random number generator for reproducible results
idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
% make simulation folder
pathSimulation = [pathCommon, 'default_network_new_seed\'];
pathResult = createresultfolder(pathSimulation);
% set new params
newParams = defaultParams;
% simulate network
S = cell(nSeed, 1);
for iSeed = 1:nSeed
    fprintf('Iteration %d: ', iSeed);
    timeStamp = tic;
    S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
    toc(timeStamp)
end % for iSeed
save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
%---------------------------------------------------------------------------------------------------%

%%%%%%%%%%--Partitioned Network Cases with Both Odors having Different Input Drives--%%%%%%%%%%
nSeed = 100;
rng('default'); % reset random number generator for reproducible results
idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
% set new params
newParams = defaultParams;
% change PN class to partitioned
newParams.networkParams.neuronLayer.vopn.class = @VectorOlfactoryProjectionNeuronLayerPartitionedTotal;
% change PN firing probability so that each odor has the same number of active PNs (required for partitioning)
newParams.simulationParams.input.vopn.fnInput = 'pseudo-random-column';
% set variation ranges
driveRange = 500:5:540;
for iSim = 1:length(driveRange)
    % add total input drive for second odor
    newParams.networkParams.neuronLayer.vopn.nTotalSpikes = driveRange(iSim);
    % make simulation folder
    pathSimulation = [pathCommon, sprintf('partitioned_network_both_odor_different_input_drive_%d\\', driveRange(iSim))];
    pathResult = createresultfolder(pathSimulation);
    % simulate network
    S = cell(nSeed, 1);
    for iSeed = 1:nSeed
        fprintf('Iteration %d: ', iSeed);
        timeStamp = tic;
        S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
        toc(timeStamp)
    end % for iSeed
    save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
end
%---------------------------------------------------------------------------------------------------%

%%%%%%%%%%--Default Network Cases with Variation in PN-KC Connection Probability--%%%%%%%%%%
nSeed = 100;
rng('default'); % reset random number generator for reproducible results
idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
% set new params
newParams = defaultParams;
% set variation ranges
cPnKcRange = [0.03 0.05 0.06 0.08 0.1:0.05:1];
for c = cPnKcRange
    % set KC number
    newParams.networkParams.layerConnection.vopn_vokc.pConnection = [c, c];
    % make simulation folder
    pathSimulation = [pathCommon, sprintf('default_network_c_pn_kc_%.2f\\', c)];
    pathResult = createresultfolder(pathSimulation);
    % simulate network
    S = cell(nSeed, 1);
    for iSeed = 1:nSeed
        fprintf('Iteration %d: ', iSeed);
        timeStamp = tic;
        S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
        toc(timeStamp)
    end % for iSeed
    save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
end
%---------------------------------------------------------------------------------------------------%

%%%%%%%%%%--Default Network Cases with Variation in KC Response Threshold--%%%%%%%%%%
nSeed = 100;
rng('default'); % reset random number generator for reproducible results
idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
% set new params
newParams = defaultParams;
% set variation ranges
thresholdRange = 50:10:250;
for thresh = thresholdRange
    % set KC number
    newParams.networkParams.neuronLayer.vokc.threshold = thresh;
    % make simulation folder
    pathSimulation = [pathCommon, sprintf('default_network_kc_thresh_%d\\', thresh)];
    pathResult = createresultfolder(pathSimulation);
    % simulate network
    S = cell(nSeed, 1);
    for iSeed = 1:nSeed
        fprintf('Iteration %d: ', iSeed);
        timeStamp = tic;
        S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
        toc(timeStamp)
    end % for iSeed
    save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
end
%---------------------------------------------------------------------------------------------------%

%%%%%%%%%%--Default Network Cases with Variation in number of PNs--%%%%%%%%%%
nSeed = 100;
rng('default'); % reset random number generator for reproducible results
idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
% set new params
newParams = defaultParams;
newParams.simulationParams.input.vopn.nInput = 100;
% set variation ranges
numRange = 20:5:100;
for nNeuron = numRange
    % set KC number
    newParams.networkParams.neuronLayer.vopn.nNeuron = nNeuron;
    % make simulation folder
    pathSimulation = [pathCommon, sprintf('default_network_num_pn_%d\\', nNeuron)];
    pathResult = createresultfolder(pathSimulation);
    % simulate network
    S = cell(nSeed, 1);
    for iSeed = 1:nSeed
        fprintf('Iteration %d: ', iSeed);
        timeStamp = tic;
        S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
        toc(timeStamp)
    end % for iSeed
    save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
end
%---------------------------------------------------------------------------------------------------%

%%%%%%%%%%--Default Network Cases with Variation in PN response probability--%%%%%%%%%%
nSeed = 100;
rng('default'); % reset random number generator for reproducible results
idSeed = sort(randperm(1e5, nSeed)); % generate seeds for the simulations
% set new params
newParams = defaultParams;
newParams.simulationParams.input.vopn.nInput = 100;
% set variation ranges
pRange = 0.1:0.05:0.9;
for p = pRange
    % set KC number
    newParams.simulationParams.input.vopn.pInput = p;
    % make simulation folder
    pathSimulation = [pathCommon, sprintf('default_network_pn_p_%.2f\\', p)];
    pathResult = createresultfolder(pathSimulation);
    % simulate network
    S = cell(nSeed, 1);
    for iSeed = 1:nSeed
        fprintf('Iteration %d: ', iSeed);
        timeStamp = tic;
        S{iSeed} = VectorOlfactorySimulation(newParams, idSeed(iSeed));
        toc(timeStamp)
    end % for iSeed
    save([pathResult, 'raw_data.mat'], 'S', 'idSeed', 'nSeed', 'defaultParams', 'newParams');
end
%---------------------------------------------------------------------------------------------------%

end % simulateallnetworkvariations