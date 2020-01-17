function runanalysisallsimulations()
%%RUNANALYSISALLSIMULATIONS Runs all analyses for all simulations in the
%results folder
%
% Usage:
%   RUNANALYSISALLSIMULATIONS()

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

pathResult = createresultfolder('mbon_stereotypy_final\');
tempNameFolder = dir(pathResult);
% remove parent folders
nameFolders = setdiff({tempNameFolder([tempNameFolder.isdir]).name}, {'.', '..', 'analyses'});
% run analyses for each folder
for eachFolder = nameFolders
    analysembonstereotypyvectorsimulation([pathResult, eachFolder{:}, '\']);
    fprintf('Analysis finished for %s\n', eachFolder{:});
end
% run special analyses
nameFolders = {'default_network_odor_100'};
for eachFolder = nameFolders
    if exist([pathResult, eachFolder{:}, '\raw_data.mat'], 'file') && ~exist([pathResult, eachFolder{:}, '\individual_stereotypy_data.mat'], 'file')
        load([pathResult, eachFolder{:}, '\raw_data.mat']);
        saveindividualstereotypyplots([pathResult, eachFolder{:}, '\'], S, nSeed);
        fprintf('Individual Analysis finished for %s\n', eachFolder{:});
    end
end
end % runanalysisallsimulations