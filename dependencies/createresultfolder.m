function resultPath = createresultfolder(nameNewFolder, overrideDefaultPath)
%%CREATERESULTFOLDER Creates and outputs the path of the simulation results
%folder appended to the current default folder
%
% Usage:
%   resultPath = CREATERESULTFOLDER(nameNewFolder, overrideDefaultPath)
%   resultPath = CREATERESULTFOLDER(nameNewFolder)
%
% Inputs:
%         nameNewFolder: name of the new results folder to be created (Outputs default result folder if no folder name specified)
%   overrideDefaultPath: specifies whether to create the new folder in the current results folder (Default: false)
%
% Output:
%   resultPath: path to the newly created folder
%
% Example:
%   CREATERESULTFOLDER('new_results')
%   Creates the 'new_results' folder (if it does not already exist) in the
%   default results folder at 'E:\aarush\neuro_lab_iitk\projects\matlab\phd\'

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

defaultPath = 'E:\aarush\neuro_lab_iitk\projects\matlab\phd\results\';
% defaultPath = [pwd, '\'];
% specify default value for the override switch
if nargin < 2
    overrideDefaultPath = false;
end
% only creates a new folder if it does not exist at the default path with
% the option of overriding the default path. Also it outputs default
% results folder if no new folder name is specified
if nargin < 1
    resultPath = defaultPath;
else
    if overrideDefaultPath
        if ~exist(nameNewFolder, 'dir')
            mkdir(nameNewFolder);
        end % if exist
        resultPath = nameNewFolder;
    else
        resultPath = [defaultPath, nameNewFolder];
        if ~exist(resultPath, 'dir')
            mkdir(resultPath);
        end % if exist
    end % if overrideDefaultPath
end % if nargin
end % createresultfolder