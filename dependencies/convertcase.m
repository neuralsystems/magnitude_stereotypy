function outputString = convertcase(inputString, varargin)
%%CONVERTCASE Convert the case of a string or a cell of strings to a
%specified case. Has capability to remove or replace all occurences of a
%specified substring
%
% Usage:
%   outputString = CONVERTCASE(inputString, outputCase, removeString, replaceString)
%   outputString = CONVERTCASE(inputString, outputCase, removeString)
%   outputString = CONVERTCASE(inputString, outputCase)
%   outputString = CONVERTCASE(inputString)
%
% Inputs:
%     inputString: string or cell of strings
%      outputCase: type of case to convert the string into (Default: 'lower')
%                  Options:
%                       'upper': UPPER CASE
%                       'lower': lower case
%                       'title': Title Case
%                    'sentence': Sentence case
%    removeString: substring to be removed (Default: '')
%   replaceString: substring which replaces the removed string (Default: '');
%
% Output:
%   outputString: string or cell of strings converted to the specified case with  specified substring removed or replaced
%
% Example:
%   CONVERTCASE('i am superman', 'title', 'superman', 'batman')
%   Outputs the string 'I Am Batman'

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

% parse inputs with specified defaults
p = inputParser();
p.addRequired('inputString', @(x) ischar(x) || iscell(x));
p.addOptional('outputCase', 'lower', @ischar);
p.addOptional('removeString', '', @ischar);
p.addOptional('replaceString', '', @ischar);
p.parse(inputString, varargin{:});
% assign parsed inputs to function variables
inputString = p.Results.inputString;
outputCase = p.Results.outputCase;
removeString = p.Results.removeString;
replaceString = p.Results.replaceString;
switch outputCase
    case 'upper'
        outputString = upper(strrep(inputString, removeString, replaceString));
    case 'lower'
        outputString = lower(strrep(inputString, removeString, replaceString));
    case 'title'
        outputString = regexprep(strrep(inputString, removeString, replaceString), '(\<[a-z])','${upper($1)}');
    case 'sentence'
        outputString = regexprep(strrep(inputString, removeString, replaceString), '(^[a-z])','${upper($1)}');
    otherwise % error for unknown case type
        disp('Unknown case type. Returning same string as before...');
        outputString = inputString;
end % switch outputCase
end % convertcase