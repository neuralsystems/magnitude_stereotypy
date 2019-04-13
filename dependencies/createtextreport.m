function createtextreport(fileData, filePath)
%%CREATETEXTREPORT Creates a text report of the specified data at the
%specified location
%
% Usage:
%   CREATETEXTREPORT(fileData, filePath)
%   CREATETEXTREPORT(fileData)
%
% Inputs:
%   fileData: string or cell of strings to be written in file
%   filePath: full path specifying the location and name of file (If unspecified, the file is written in the current active folder with the filename 'new_report.txt'

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

% write a new file ensuring the default filePath
if nargin < 2
    filePath = [pwd, '\new_report.txt'];
end
fid = fopen(filePath, 'w');
if ~(fid == 1)
        fprintf(fid, '%s', fileData{:});
end
fclose(fid);
end % createtextreport