function computevaluesschaffercode()
%%COMPUTEVALUESSCHAFFERCODE Calculates the stereotypy and correlation
%values for the untrained readout neuron from the schaffer code results.
%Need to run the calculateSumFigParts.m script 6 times before running this
%
% Usage:
%   COMPUTEVALUESSCHAFFERCODE()

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

nSim = 6;
dataStereotypyUnmod = zeros(nSim, 3);
dataStereotypyUnmodTrained = zeros(nSim, 3);
dataStereotypyNoMean = zeros(nSim, 3);
dataStereotypyNoMeanTrained = zeros(nSim, 3);
dataCorrUnmod = zeros(nSim, 3);
dataCorrUnmodTrained = zeros(nSim, 3);
dataCorrNoMean = zeros(nSim, 3);
dataCorrNoMeanTrained = zeros(nSim, 3);
for iSim = 1:nSim
    %%%% correlation and stereotypy in readouts of unmodified code
    load(['data/run_', num2str(iSim), '/partsForSumFig.mat'], 'UNTRAINEDFourthOrder90p*', 'trainedFourthOrder90p*', 'trainedOdorNum')
    normQt1 = trainedFourthOrder90pClass1(trainedOdorNum);
    normQt2 = trainedFourthOrder90pClass1_mouse2(trainedOdorNum);
    
    %%%% Untrained Data %%%%
    % Class1 (f = 0.7)
    data.class1 = full([UNTRAINEDFourthOrder90pClass1 / normQt1; UNTRAINEDFourthOrder90pClass1_mouse2 / normQt2]);
    stereotypy.class1 = computeindividualpairstereotypy(data.class1);
    corrCoef.class1 = corrcoef(data.class1.');
    % Class2 (f = 0.3)
    data.class2 = full([UNTRAINEDFourthOrder90pClass2 / normQt1; UNTRAINEDFourthOrder90pClass2_mouse2 / normQt2]);
    stereotypy.class2 = computeindividualpairstereotypy(data.class2);
    corrCoef.class2 = corrcoef(data.class2.');
    % Class3 (f = 0)
    data.class3 = full([UNTRAINEDFourthOrder90pRand / normQt1; UNTRAINEDFourthOrder90pRand_mouse2 / normQt2]);
    stereotypy.class3 = computeindividualpairstereotypy(data.class3);
    corrCoef.class3 = corrcoef(data.class3.');
    % store data
    dataStereotypyUnmod(iSim, :) = [mean(stereotypy.class1), mean(stereotypy.class2), mean(stereotypy.class3)];
    dataCorrUnmod(iSim, :) = [corrCoef.class1(2), corrCoef.class2(2), corrCoef.class3(2)];
    
    %%%% Trained Data %%%%
    % Class1 (f = 0.7)
    data.class1 = full([trainedFourthOrder90pClass1 / normQt1; trainedFourthOrder90pClass1_mouse2 / normQt2]);
    stereotypy.class1 = computeindividualpairstereotypy(data.class1);
    corrCoef.class1 = corrcoef(data.class1.');
    % Class2 (f = 0.3)
    data.class2 = full([trainedFourthOrder90pClass2 / normQt1; trainedFourthOrder90pClass2_mouse2 / normQt2]);
    stereotypy.class2 = computeindividualpairstereotypy(data.class2);
    corrCoef.class2 = corrcoef(data.class2.');
    % Class3 (f = 0)
    data.class3 = full([trainedFourthOrder90pRand / normQt1; trainedFourthOrder90pRand_mouse2 / normQt2]);
    stereotypy.class3 = computeindividualpairstereotypy(data.class3);
    corrCoef.class3 = corrcoef(data.class3.');
    % store data
    dataStereotypyUnmodTrained(iSim, :) = [mean(stereotypy.class1), mean(stereotypy.class2), mean(stereotypy.class3)];
    dataCorrUnmodTrained(iSim, :) = [corrCoef.class1(2), corrCoef.class2(2), corrCoef.class3(2)];
    
    %%%% correlation in untrained readouts of modified code
    load(['data/run_', num2str(iSim), '/partsForSumFigNoMean.mat'], 'UNTRAINEDFourthOrder90p*', 'trainedFourthOrder90p*', 'trainedOdorNum')
    normQt1 = trainedFourthOrder90pClass1(trainedOdorNum);
    normQt2 = trainedFourthOrder90pClass1_mouse2(trainedOdorNum);
    
    %%%% Untrained Data %%%%
    % Class1 (f = 0.7)
    data.class1 = full([UNTRAINEDFourthOrder90pClass1 / normQt1; UNTRAINEDFourthOrder90pClass1_mouse2 / normQt2]);
    stereotypy.class1 = computeindividualpairstereotypy(data.class1);
    corrCoef.class1 = corrcoef(data.class1.');
    % Class2 (f = 0.3)
    data.class2 = full([UNTRAINEDFourthOrder90pClass2 / normQt1; UNTRAINEDFourthOrder90pClass2_mouse2 / normQt2]);
    stereotypy.class2 = computeindividualpairstereotypy(data.class2);
    corrCoef.class2 = corrcoef(data.class2.');
    % Class3 (f = 0)
    data.class3 = full([UNTRAINEDFourthOrder90pRand / normQt1; UNTRAINEDFourthOrder90pRand_mouse2 / normQt2]);
    stereotypy.class3 = computeindividualpairstereotypy(data.class3);
    corrCoef.class3 = corrcoef(data.class3.');
    % store data
    dataStereotypyNoMean(iSim, :) = [mean(stereotypy.class1), mean(stereotypy.class2), mean(stereotypy.class3)];
    dataCorrNoMean(iSim, :) = [corrCoef.class1(2), corrCoef.class2(2), corrCoef.class3(2)];
    
    %%%% Trained Data %%%%
    % Class1 (f = 0.7)
    data.class1 = full([trainedFourthOrder90pClass1 / normQt1; trainedFourthOrder90pClass1_mouse2 / normQt2]);
    stereotypy.class1 = computeindividualpairstereotypy(data.class1);
    corrCoef.class1 = corrcoef(data.class1.');
    % Class2 (f = 0.3)
    data.class2 = full([trainedFourthOrder90pClass2 / normQt1; trainedFourthOrder90pClass2_mouse2 / normQt2]);
    stereotypy.class2 = computeindividualpairstereotypy(data.class2);
    corrCoef.class2 = corrcoef(data.class2.');
    % Class3 (f = 0)
    data.class3 = full([trainedFourthOrder90pRand / normQt1; trainedFourthOrder90pRand_mouse2 / normQt2]);
    stereotypy.class3 = computeindividualpairstereotypy(data.class3);
    corrCoef.class3 = corrcoef(data.class3.');
    % store data
    dataStereotypyNoMeanTrained(iSim, :) = [mean(stereotypy.class1), mean(stereotypy.class2), mean(stereotypy.class3)];
    dataCorrNoMeanTrained(iSim, :) = [corrCoef.class1(2), corrCoef.class2(2), corrCoef.class3(2)];
    disp(iSim)
end % for iSim
save('data/schaffer_code_data.mat', 'dataStereotypyUnmod', 'dataStereotypyUnmodTrained','dataStereotypyNoMean', 'dataStereotypyNoMeanTrained', 'dataCorrUnmod', 'dataCorrUnmodTrained', 'dataCorrNoMean', 'dataCorrNoMeanTrained');
end % computevaluesschaffercode