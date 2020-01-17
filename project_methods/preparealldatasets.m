function preparealldatasets()

% Stereotypy data for theoretical analysis
if ~exist('data\theoretical_stereotypy\theoretical_stereotypy.mat', 'file')
    computetheoreticaltotalkcresponsestereotypy();
end
fprintf('Processed theoretical stereotypy data...\n')

% Stereotypy data for drosophila gcamp3 recordings
raw_folder = 'data\fly_gcamp3\';
processed_file = 'data\fly_gcamp3\fly_data.mat';
if ~exist(processed_file, 'file')
    % Import the data
    [~, ~, raw] = xlsread([raw_folder 'fly_mbon_gcamp3_stereotypy.xlsx'], 'stereotypies - 3 odors', 'A5:AD634');
    raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x), raw)) = {''};
    % Replace non-numeric cells with NaN
    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x), raw); % Find non-numeric cells
    raw(R) = {NaN}; % Replace non-numeric cells
    % Create output variable
    data = reshape([raw{:}], size(raw));
    % Allocate imported array to variable names for each lobe
    control.alpha = reshape(data(:, 13:15), [], 1);
    control.alpha(isnan(control.alpha)) = [];
    apl_tnt.alpha = reshape(data(:, 16:18), [], 1);
    apl_tnt.alpha(isnan(apl_tnt.alpha)) = [];
    control.alpha_p = reshape(data(:, 1:3), [], 1);
    control.alpha_p(isnan(control.alpha_p)) = [];
    apl_tnt.alpha_p = reshape(data(:, 4:6), [], 1);
    apl_tnt.alpha_p(isnan(apl_tnt.alpha_p)) = [];
    control.beta = reshape(data(:, 19:21), [], 1);
    control.beta(isnan(control.beta)) = [];
    apl_tnt.beta = reshape(data(:, 22:24), [], 1);
    apl_tnt.beta(isnan(apl_tnt.beta)) = [];
    control.beta_p = reshape(data(:, 7:9), [], 1);
    control.beta_p(isnan(control.beta_p)) = [];
    apl_tnt.beta_p = reshape(data(:, 10:12), [], 1);
    apl_tnt.beta_p(isnan(apl_tnt.beta_p)) = [];
    control.gamma = reshape(data(:, 25:27), [], 1);
    control.gamma(isnan(control.gamma)) = [];
    apl_tnt.gamma = reshape(data(:, 28:30), [], 1);
    apl_tnt.gamma(isnan(apl_tnt.gamma)) = [];
    save(processed_file, 'control', 'apl_tnt')
    clearvars
end
fprintf('Processed fly gcamp3 data...\n')

% Schaffer, E. S. et al. Odor Perception on the Two Sides of the Brain: Consistency Despite Randomness. Neuron 98, 736-742.e3 (2018).
raw_folder = 'data\schaffer_2018\raw\';
stereo_file = 'data\schaffer_2018\stereotypy_data_schaffer.mat';
if ~exist(stereo_file, 'file')
    nRun = 6;
    dataCorrNoMean = zeros(nRun, 3);
    dataCorrNoMeanTrained = zeros(nRun, 3);
    dataCorrUnmod = zeros(nRun, 3);
    dataCorrUnmodTrained = zeros(nRun, 3);
    dataStereotypyNoMean = zeros(nRun, 3);
    dataStereotypyNoMeanTrained = zeros(nRun, 3);
    dataStereotypyUnmod = zeros(nRun, 3);
    dataStereotypyUnmodTrained = zeros(nRun, 3);
    nClass = 1000;
    nRand = 2000;
    idCombinationsClass = nchoosek(1:nClass, 2);
    idCombinationsRand = nchoosek(1:nRand, 2);
    for iRun = 1:nRun
        % data with weight normalization
        load(sprintf('%srun_%d\\partsForSumFig.mat', raw_folder, iRun), 'UNTRAINEDFourthOrder90p*', 'trainedFourthOrder90p*')
        % Class1 (f = 0.7)
        [dataStereotypyUnmodTrained(iRun, 1), dataCorrUnmodTrained(iRun, 1)] = computeschaffercase(trainedFourthOrder90pClass1, trainedFourthOrder90pClass1_mouse2, idCombinationsClass);
        [dataStereotypyUnmod(iRun, 1), dataCorrUnmod(iRun, 1)] = computeschaffercase(UNTRAINEDFourthOrder90pClass1, UNTRAINEDFourthOrder90pClass1_mouse2, idCombinationsClass);
        % Class2 (f = 0.3)
        [dataStereotypyUnmodTrained(iRun, 2), dataCorrUnmodTrained(iRun, 2)] = computeschaffercase(trainedFourthOrder90pClass2, trainedFourthOrder90pClass2_mouse2, idCombinationsClass);
        [dataStereotypyUnmod(iRun, 2), dataCorrUnmod(iRun, 2)] = computeschaffercase(UNTRAINEDFourthOrder90pClass2, UNTRAINEDFourthOrder90pClass2_mouse2, idCombinationsClass);
        % Class3 (f = 0)
        [dataStereotypyUnmodTrained(iRun, 3), dataCorrUnmodTrained(iRun, 3)] = computeschaffercase(trainedFourthOrder90pRand, trainedFourthOrder90pRand_mouse2, idCombinationsRand);
        [dataStereotypyUnmod(iRun, 3), dataCorrUnmod(iRun, 3)] = computeschaffercase(UNTRAINEDFourthOrder90pRand, UNTRAINEDFourthOrder90pRand_mouse2, idCombinationsRand);
        % data without weight normalization
        load(sprintf('%srun_%d\\partsForSumFigNoMean.mat', raw_folder, iRun), 'UNTRAINEDFourthOrder90p*', 'trainedFourthOrder90p*')
        % Class1 (f = 0.7)
        [dataStereotypyNoMeanTrained(iRun, 1), dataCorrNoMeanTrained(iRun, 1)] = computeschaffercase(trainedFourthOrder90pClass1, trainedFourthOrder90pClass1_mouse2, idCombinationsClass);
        [dataStereotypyNoMean(iRun, 1), dataCorrNoMean(iRun, 1)] = computeschaffercase(UNTRAINEDFourthOrder90pClass1, UNTRAINEDFourthOrder90pClass1_mouse2, idCombinationsClass);
        % Class2 (f = 0.3)
        [dataStereotypyNoMeanTrained(iRun, 2), dataCorrNoMeanTrained(iRun, 2)] = computeschaffercase(trainedFourthOrder90pClass2, trainedFourthOrder90pClass2_mouse2, idCombinationsClass);
        [dataStereotypyNoMean(iRun, 2), dataCorrNoMean(iRun, 2)] = computeschaffercase(UNTRAINEDFourthOrder90pClass2, UNTRAINEDFourthOrder90pClass2_mouse2, idCombinationsClass);
        % Class3 (f = 0)
        [dataStereotypyNoMeanTrained(iRun, 3), dataCorrNoMeanTrained(iRun, 3)] = computeschaffercase(trainedFourthOrder90pRand, trainedFourthOrder90pRand_mouse2, idCombinationsRand);
        [dataStereotypyNoMean(iRun, 3), dataCorrNoMean(iRun, 3)] = computeschaffercase(UNTRAINEDFourthOrder90pRand, UNTRAINEDFourthOrder90pRand_mouse2, idCombinationsRand);
    end
    save(stereo_file, 'data*');
    clearvars
end
fprintf('Processed schaffer et al. 2018 data...\n')

% Shimizu, K. & Stopfer, M. A Population of Projection Neurons that Inhibits the Lateral Horn but Excites the Antennal Lobe through Chemical Synapses in Drosophila. Front. Neural Circuits 11, (2017).
raw_folder = 'data\shimizu_2017\raw\';
processed_file = 'data\shimizu_2017\processed_data_shimizu.mat';
if ~exist(processed_file, 'file')
    info_data.id_odor = {'benzaldehyde', '2-octanone', 'ethyl acetate', 'pentylacetate', 'ethyl butyrate'};
    info_data.id_glomerulus = {'VC4', 'DL2v', 'VM5v', 'VC3'};
    
    % Convert data to nPn (glomerulus) x nTrial x nIndividual (brain) x nOdor. Sum over the time frames
    processed_data = nan(4, 10, 6, 5); % size decided using maximum observed values
    response_time_start = 2;
    response_time_end = 4;
    
    % vc4 data
    iPn = 1;
    relevant_odors = 1:5;
    
    load([raw_folder '130211.mat'])
    current_data = spikeTMatrix(1:35, :);
    temp_data = (sum(current_data >= response_time_start & current_data <= response_time_end, 2) / (response_time_end - response_time_start)) - (sum(current_data < response_time_start & current_data ~= 0, 2) / response_time_start);
    odor_start_indices = [1 8 15 22 29];
    odor_end_indices = [7 14 21 28 35];
    for iOdor = 1:length(odor_start_indices)
        processed_data(iPn, 1:(odor_end_indices(iOdor) - odor_start_indices(iOdor) + 1), 1, relevant_odors(iOdor)) = temp_data(odor_start_indices(iOdor):odor_end_indices(iOdor), 1);
    end
    
    load([raw_folder '130214.mat'])
    current_data = spikeTMatrix(1:32, :);
    temp_data = (sum(current_data >= response_time_start & current_data <= response_time_end, 2) / (response_time_end - response_time_start)) - (sum(current_data < response_time_start & current_data ~= 0, 2) / response_time_start);
    odor_start_indices = [1 7 13 19 26];
    odor_end_indices = [6 12 18 25 32];
    for iOdor = 1:length(odor_start_indices)
        processed_data(iPn, 1:(odor_end_indices(iOdor) - odor_start_indices(iOdor) + 1), 2, relevant_odors(iOdor)) = temp_data(odor_start_indices(iOdor):odor_end_indices(iOdor), 1);
    end
    
    % dl2v data
    iPn = 2;
    relevant_odors = [1 5];
    
    load([raw_folder '151225.mat'])
    current_data = spikeTMatrix(1:17, :);
    temp_data = (sum(current_data >= response_time_start & current_data <= response_time_end, 2) / (response_time_end - response_time_start)) - (sum(current_data < response_time_start & current_data ~= 0, 2) / response_time_start);
    odor_start_indices = [1 8];
    odor_end_indices = [7 17];
    for iOdor = 1:length(odor_start_indices)
        processed_data(iPn, 1:(odor_end_indices(iOdor) - odor_start_indices(iOdor) + 1), 1, relevant_odors(iOdor)) = temp_data(odor_start_indices(iOdor):odor_end_indices(iOdor), 1);
    end
    
    load([raw_folder '160101.mat'])
    current_data = spikeTMatrix(1:21, :);
    temp_data = (sum(current_data >= response_time_start & current_data <= response_time_end, 2) / (response_time_end - response_time_start)) - (sum(current_data < response_time_start & current_data ~= 0, 2) / response_time_start);
    odor_start_indices = [1 15];
    odor_end_indices = [7 21];
    for iOdor = 1:length(odor_start_indices)
        processed_data(iPn, 1:(odor_end_indices(iOdor) - odor_start_indices(iOdor) + 1), 2, relevant_odors(iOdor)) = temp_data(odor_start_indices(iOdor):odor_end_indices(iOdor), 1);
    end
    
    load([raw_folder '160105-1.mat'])
    current_data = spikeTMatrix(1:21, :);
    temp_data = (sum(current_data >= response_time_start & current_data <= response_time_end, 2) / (response_time_end - response_time_start)) - (sum(current_data < response_time_start & current_data ~= 0, 2) / response_time_start);
    odor_start_indices = [1 15];
    odor_end_indices = [7 21];
    for iOdor = 1:length(odor_start_indices)
        processed_data(iPn, 1:(odor_end_indices(iOdor) - odor_start_indices(iOdor) + 1), 3, relevant_odors(iOdor)) = temp_data(odor_start_indices(iOdor):odor_end_indices(iOdor), 1);
    end
    
    load([raw_folder '151101.mat'])
    current_data = spikeTMatrix(1:14, :);
    temp_data = (sum(current_data >= response_time_start & current_data <= response_time_end, 2) / (response_time_end - response_time_start)) - (sum(current_data < response_time_start & current_data ~= 0, 2) / response_time_start);
    odor_start_indices = [1 8];
    odor_end_indices = [7 14];
    for iOdor = 1:length(odor_start_indices)
        processed_data(iPn, 1:(odor_end_indices(iOdor) - odor_start_indices(iOdor) + 1), 4, relevant_odors(iOdor)) = temp_data(odor_start_indices(iOdor):odor_end_indices(iOdor), 1);
    end
    
    load([raw_folder '160116-2.mat'])
    current_data = spikeTMatrix(1:28, :);
    temp_data = (sum(current_data >= response_time_start & current_data <= response_time_end, 2) / (response_time_end - response_time_start)) - (sum(current_data < response_time_start & current_data ~= 0, 2) / response_time_start);
    odor_start_indices = [1 22];
    odor_end_indices = [7 28];
    for iOdor = 1:length(odor_start_indices)
        processed_data(iPn, 1:(odor_end_indices(iOdor) - odor_start_indices(iOdor) + 1), 5, relevant_odors(iOdor)) = temp_data(odor_start_indices(iOdor):odor_end_indices(iOdor), 1);
    end
    
    load([raw_folder '160118.mat'])
    current_data = spikeTMatrix(1:7, :);
    temp_data = (sum(current_data >= response_time_start & current_data <= response_time_end, 2) / (response_time_end - response_time_start)) - (sum(current_data < response_time_start & current_data ~= 0, 2) / response_time_start);
    odor_start_indices = 1;
    odor_end_indices = 7;
    for iOdor = 1:length(odor_start_indices)
        processed_data(iPn, 1:(odor_end_indices(iOdor) - odor_start_indices(iOdor) + 1), 6, relevant_odors(iOdor)) = temp_data(odor_start_indices(iOdor):odor_end_indices(iOdor), 1);
    end
    % odor 2 extracted manually
    spike_counts_response = [88 75 72 72 68 65 71];
    spike_counts_background = [7 0 0 2 0 6 2];
    processed_data(iPn, 1:length(spike_counts_response), 6, relevant_odors(2)) = (spike_counts_response / (response_time_end - response_time_start)) - (spike_counts_background / response_time_start);
    
    % vm5v data
    iPn = 3;
    relevant_odors = [2 4 5];
    
    load([raw_folder '151226-1.mat'])
    current_data = spikeTMatrix(1:21, :);
    temp_data = (sum(current_data >= response_time_start & current_data <= response_time_end, 2) / (response_time_end - response_time_start)) - (sum(current_data < response_time_start & current_data ~= 0, 2) / response_time_start);
    odor_start_indices = [1 8 15];
    odor_end_indices = [7 14 21];
    for iOdor = 1:length(odor_start_indices)
        processed_data(iPn, 1:(odor_end_indices(iOdor) - odor_start_indices(iOdor) + 1), 1, relevant_odors(iOdor)) = temp_data(odor_start_indices(iOdor):odor_end_indices(iOdor), 1);
    end
    
    load([raw_folder '151229.mat'])
    current_data = spikeTMatrix(1:21, :);
    temp_data = (sum(current_data >= response_time_start & current_data <= response_time_end, 2) / (response_time_end - response_time_start)) - (sum(current_data < response_time_start & current_data ~= 0, 2) / response_time_start);
    odor_start_indices = [1 8 15];
    odor_end_indices = [7 14 21];
    for iOdor = 1:length(odor_start_indices)
        processed_data(iPn, 1:(odor_end_indices(iOdor) - odor_start_indices(iOdor) + 1), 2, relevant_odors(iOdor)) = temp_data(odor_start_indices(iOdor):odor_end_indices(iOdor), 1);
    end
    
    load([raw_folder '160104-2.mat'])
    current_data = spikeTMatrix(1:21, :);
    temp_data = (sum(current_data >= response_time_start & current_data <= response_time_end, 2) / (response_time_end - response_time_start)) - (sum(current_data < response_time_start & current_data ~= 0, 2) / response_time_start);
    odor_start_indices = [1 8 15];
    odor_end_indices = [7 14 21];
    for iOdor = 1:length(odor_start_indices)
        processed_data(iPn, 1:(odor_end_indices(iOdor) - odor_start_indices(iOdor) + 1), 3, relevant_odors(iOdor)) = temp_data(odor_start_indices(iOdor):odor_end_indices(iOdor), 1);
    end
    
    load([raw_folder '160105-2.mat'])
    current_data = spikeTMatrix(1:21, :);
    temp_data = (sum(current_data >= response_time_start & current_data <= response_time_end, 2) / (response_time_end - response_time_start)) - (sum(current_data < response_time_start & current_data ~= 0, 2) / response_time_start);
    odor_start_indices = [1 8 15];
    odor_end_indices = [7 14 21];
    for iOdor = 1:length(odor_start_indices)
        processed_data(iPn, 1:(odor_end_indices(iOdor) - odor_start_indices(iOdor) + 1), 4, relevant_odors(iOdor)) = temp_data(odor_start_indices(iOdor):odor_end_indices(iOdor), 1);
    end
    
    load([raw_folder '160331.mat'])
    current_data = spikeTMatrix(1:21, :);
    temp_data = (sum(current_data >= response_time_start & current_data <= response_time_end, 2) / (response_time_end - response_time_start)) - (sum(current_data < response_time_start & current_data ~= 0, 2) / response_time_start);
    odor_start_indices = [1 8 15];
    odor_end_indices = [7 14 21];
    for iOdor = 1:length(odor_start_indices)
        processed_data(iPn, 1:(odor_end_indices(iOdor) - odor_start_indices(iOdor) + 1), 5, relevant_odors(iOdor)) = temp_data(odor_start_indices(iOdor):odor_end_indices(iOdor), 1);
    end
    
    load([raw_folder '160404.mat'])
    current_data = spikeTMatrix(1:21, :);
    temp_data = (sum(current_data >= response_time_start & current_data <= response_time_end, 2) / (response_time_end - response_time_start)) - (sum(current_data < response_time_start & current_data ~= 0, 2) / response_time_start);
    odor_start_indices = [1 8 15];
    odor_end_indices = [7 14 21];
    for iOdor = 1:length(odor_start_indices)
        processed_data(iPn, 1:(odor_end_indices(iOdor) - odor_start_indices(iOdor) + 1), 6, relevant_odors(iOdor)) = temp_data(odor_start_indices(iOdor):odor_end_indices(iOdor), 1);
    end
    
    % vc3 data
    iPn = 4;
    relevant_odors = [1 4 5];
    
    load([raw_folder '160329.mat'])
    current_data = spikeTMatrix(1:28, :);
    temp_data = (sum(current_data >= response_time_start & current_data <= response_time_end, 2) / (response_time_end - response_time_start)) - (sum(current_data < response_time_start & current_data ~= 0, 2) / response_time_start);
    odor_start_indices = [22 1 8];
    odor_end_indices = [28 7 14];
    for iOdor = 1:length(odor_start_indices)
        processed_data(iPn, 1:(odor_end_indices(iOdor) - odor_start_indices(iOdor) + 1), 1, relevant_odors(iOdor)) = temp_data(odor_start_indices(iOdor):odor_end_indices(iOdor), 1);
    end
    
    load([raw_folder '160506-2.mat'])
    current_data = spikeTMatrix(1:21, :);
    temp_data = (sum(current_data >= response_time_start & current_data <= response_time_end, 2) / (response_time_end - response_time_start)) - (sum(current_data < response_time_start & current_data ~= 0, 2) / response_time_start);
    odor_start_indices = [1 8 15];
    odor_end_indices = [7 14 21];
    for iOdor = 1:length(odor_start_indices)
        processed_data(iPn, 1:(odor_end_indices(iOdor) - odor_start_indices(iOdor) + 1), 2, relevant_odors(iOdor)) = temp_data(odor_start_indices(iOdor):odor_end_indices(iOdor), 1);
    end
    
    load([raw_folder '160420-2.mat'])
    current_data = spikeTMatrix(1:21, :);
    temp_data = (sum(current_data >= response_time_start & current_data <= response_time_end, 2) / (response_time_end - response_time_start)) - (sum(current_data < response_time_start & current_data ~= 0, 2) / response_time_start);
    odor_start_indices = [8 1 15];
    odor_end_indices = [14 7 21];
    for iOdor = 1:length(odor_start_indices)
        processed_data(iPn, 1:(odor_end_indices(iOdor) - odor_start_indices(iOdor) + 1), 3, relevant_odors(iOdor)) = temp_data(odor_start_indices(iOdor):odor_end_indices(iOdor), 1);
    end
    
    load([raw_folder '160430.mat'])
    current_data = spikeTMatrix(1:21, :);
    temp_data = (sum(current_data >= response_time_start & current_data <= response_time_end, 2) / (response_time_end - response_time_start)) - (sum(current_data < response_time_start & current_data ~= 0, 2) / response_time_start);
    odor_start_indices = [8 1 15];
    odor_end_indices = [14 7 21];
    for iOdor = 1:length(odor_start_indices)
        processed_data(iPn, 1:(odor_end_indices(iOdor) - odor_start_indices(iOdor) + 1), 4, relevant_odors(iOdor)) = temp_data(odor_start_indices(iOdor):odor_end_indices(iOdor), 1);
    end
    
    load([raw_folder '160116-1.mat'])
    current_data = spikeTMatrix(1:7, :);
    temp_data = (sum(current_data >= response_time_start & current_data <= response_time_end, 2) / (response_time_end - response_time_start)) - (sum(current_data < response_time_start & current_data ~= 0, 2) / response_time_start);
    odor_start_indices = 1;
    odor_end_indices = 7;
    for iOdor = 1:length(odor_start_indices)
        processed_data(iPn, 1:(odor_end_indices(iOdor) - odor_start_indices(iOdor) + 1), 5, relevant_odors(iOdor)) = temp_data(odor_start_indices(iOdor):odor_end_indices(iOdor), 1);
    end
    % odor 2 extracted manually
    spike_counts_response = [226 215 210 215 215 215 218];
    spike_counts_background = [7 7 8 1 5 6 5];
    processed_data(iPn, 1:length(spike_counts_response), 5, relevant_odors(2)) = (spike_counts_response / (response_time_end - response_time_start)) - (spike_counts_background / response_time_start);
    % odor 3 extracted manually
    spike_counts_response = [55 56 59 60 58 65 67];
    spike_counts_background = [11 20 26 22 21 12 16];
    processed_data(iPn, 1:length(spike_counts_response), 5, relevant_odors(3)) = (spike_counts_response / (response_time_end - response_time_start)) - (spike_counts_background / response_time_start);
    
    load([raw_folder '160119-1.mat'])
    current_data = spikeTMatrix(1:21, :);
    temp_data = (sum(current_data >= response_time_start & current_data <= response_time_end, 2) / (response_time_end - response_time_start)) - (sum(current_data < response_time_start & current_data ~= 0, 2) / response_time_start);
    odor_start_indices = [1 8 15];
    odor_end_indices = [7 14 21];
    for iOdor = 1:length(odor_start_indices)
        processed_data(iPn, 1:(odor_end_indices(iOdor) - odor_start_indices(iOdor) + 1), 6, relevant_odors(iOdor)) = temp_data(odor_start_indices(iOdor):odor_end_indices(iOdor), 1);
    end
    save(processed_file, 'info_data', 'processed_data');
    stereo_file = 'data\shimizu_2017\stereotypy_data_shimizu.mat';
    if ~exist(stereo_file, 'file')
        nPn = size(processed_data, 1);
        stereo_io_pred_raw = cell(nPn, 1);
        stereo_io_corr_raw = cell(nPn, 1);
        for iPn = 1:nPn
            % Stereotypy IO
            temp_data = removeallnandims(squeeze(nanmean(processed_data(iPn, :, :, :), 2)));
            if isvector(temp_data) % ignore data where only 1 individual or 1 odor exists
                continue
            end
            stereo_io_pred_raw{iPn} = stripnans(computeindividualpairstereotypy(temp_data, false));
            stereo_io_corr_raw{iPn} = stripnans(computeindividualpairpearsoncorrelation(temp_data));
        end
        stereo_io_pred_mean = cellfun(@nanmean, stereo_io_pred_raw);
        stereo_io_corr_mean = cellfun(@nanmean, stereo_io_corr_raw);
        save(stereo_file, 'stereo*');
    end
    clearvars
end
fprintf('Processed shimizu & stopfer 2017 data...\n')

% Murthy, M., Fiete, I. & Laurent, G. Testing Odor Response Stereotypy in the Drosophila Mushroom Body. Neuron 59, 1009–1023 (2008).
raw_folder = 'data\murthy_2008\';
% all odor data
processed_file = 'data\murthy_2008\processed_data_murthy.mat';
if ~exist(processed_file, 'file')
    [~, ~, raw_data] = xlsread([raw_folder 'murthy_2008.xlsx'],'Sheet1','B3:AM52');
    raw_data(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x), raw_data)) = {''};
    raw_data(cellfun(@(x) ~isnumeric(x) && ~islogical(x), raw_data)) = {NaN}; % Replace non-numeric cells
    
    % Create output variable
    raw_data = reshape([raw_data{:}],size(raw_data));
    [~, ~, info_data.id_odor] = xlsread([raw_folder 'murthy_2008.xlsx'],'Sheet1','B2:AM2');
    info_data.id_odor(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x), info_data.id_odor)) = {''};
    
    [~, ~, info_data.id_glomerulus] = xlsread([raw_folder 'murthy_2008.xlsx'],'Sheet1','A3:A52');
    info_data.id_glomerulus(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x), info_data.id_glomerulus)) = {''};
    [info_data.id_glomerulus, ~, count] = unique(info_data.id_glomerulus, 'stable');
    info_data.id_glomerulus = info_data.id_glomerulus.';
    count = histcounts(count);
    start_ind = cumsum([1 count(1:end-1)]);
    end_ind = cumsum(count);
    % Convert data to nKc (class) x nTrial x nIndividual (brain) x nOdor. Sum over the time frames
    nIndividual = max(count);
    nKc = length(count);
    nOdor = length(info_data.id_odor);
    nTrial = 1;
    processed_data = nan(nKc, nTrial, nIndividual, nOdor); % size decided using maximum observed values
    for iKc = 1:nKc
        processed_data(iKc, :, 1:count(iKc), :) = raw_data(start_ind(iKc):end_ind(iKc), :);
    end
    save(processed_file, 'info_data', 'processed_data');
    stereo_file = 'data\murthy_2008\stereotypy_data_murthy.mat';
    if ~exist(stereo_file, 'file')
        nKc = size(processed_data, 1);
        stereo_io_pred_raw = cell(nKc + 1, 1);
        stereo_io_corr_raw = cell(nKc + 1, 1);
        for iKc = 1:nKc
            % Stereotypy IO
            temp_data = removeallnandims(squeeze(nanmean(processed_data(iKc, :, :, :), 2)));
            if isvector(temp_data) % ignore data where only 1 individual or 1 odor exists
                continue
            end
            stereo_io_pred_raw{iKc} = stripnans(computeindividualpairstereotypy(temp_data, false));
            stereo_io_corr_raw{iKc} = stripnans(computeindividualpairpearsoncorrelation(temp_data));
        end
        temp_data = removeallnandims(cell2mat(arrayfun(@(x) squeeze(processed_data(x, 1, :, :)), 1:nKc, 'UniformOutput', 0).'));
        stereo_io_pred_raw{end} = stripnans(computeindividualpairstereotypy(temp_data, false));
        stereo_io_corr_raw{end} = stripnans(computeindividualpairpearsoncorrelation(temp_data));
        stereo_io_pred_mean = cellfun(@nanmean, stereo_io_pred_raw);
        stereo_io_corr_mean = cellfun(@nanmean, stereo_io_corr_raw);
        save(stereo_file, 'stereo*');
    end
    clearvars
end
fprintf('Processed murthy et al. 2008 data...\n')
end

function data = removeallnandims(data)
nDims = ndims(data);
    function allnandims = findalldim(data, dims)
        if isempty(dims)
            allnandims = data;
        else
            allnandims = findalldim(all(data, dims(1)), dims(2:end));
        end
    end
for iDim = 1:nDims
    subs_data = repelem({':'}, 1, nDims);
    subs_data{iDim} = ~findalldim(isnan(data), setdiff(1:nDims, iDim));
    data = subsref(data, substruct('()', subs_data));
end
end

function data = stripnans(data)
data(isnan(data)) = [];
data = reshape(data, [], 1);
end

function [stereo, pcorr] = computeschaffercase(ind1, ind2, idComb)
nComb = size(idComb, 1);
stereotypy = zeros(nComb, 1);
for iComb = 1:nComb
    d1 = (ind1(idComb(iComb, 1)) - ind2(idComb(iComb, 1))) ^ 2 + (ind1(idComb(iComb, 2)) - ind2(idComb(iComb, 2))) ^ 2;
    d2 = (ind1(idComb(iComb, 2)) - ind2(idComb(iComb, 1))) ^ 2 + (ind1(idComb(iComb, 1)) - ind2(idComb(iComb, 2))) ^ 2;
    stereotypy(iComb) = (d2 - d1) / (d2 + d1);
end
pcorr = corr(ind1.', ind2.');
stereo = mean(stereotypy);
end