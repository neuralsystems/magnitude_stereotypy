function main()
%%MAIN Generates all the simulations and analyses for the project. Make
%sure that the mbon_stereotypy_project folder is the current MATLAB folder
%
% Usage:
%   MAIN()

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

%%%% run all simulations
simulateallnetworkvariations

%%%% analyse all individual simulations
runanalysisallsimulations

%%%% generate simulation figures
generatefigures

end % main