function computetheoreticaltotalkcresponsestereotypy()
%%COMPUTETHEORETICALTOTALKCRESPONSESTEREOTYPY Computes the expected stereotypy
%between the sum of KC reponses for different odors and individual
%combinations. Here the PN response is assumed to be a binary vector of
%size nP with Pr[1] = p and the connection matrix to be a binary matrix of
%size nK x nP with Pr[1] = c. The distance is found over different
%conditions on PN response vectors. The KC responses are also binary
%
% !!!REQUIRES LARGE AMOUNT OF RAM TO RUN (AT LEAST 250 GB)!!!
%
% Usage:
%   COMPUTETHEORETICALTOTALKCRESPONSESTEREOTYPY()

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

nP = 50; % number of PNs
c = 0.14; % probability of 1 in connection matrix
p = 0.5; % probability of 1 in response matrix
q = 0.1; % percentage of active KCs
t = find(binocdf(0:nP, nP, p * c) <= (1 - q), 1, 'last') + 1; % threshold
nKRange = [100 250 500 1000:1000:9000]; % Range of nK values

% define stereotypy
fnStereotypy = @(x, y) (x - y) ./ (x + y);

% variation with nK
D1 = zeros(1, length(nKRange));
D2 = zeros(1, length(nKRange));
for indexM = 1:length(nKRange)
    % calculate distances
    [D1(indexM), D2(indexM)] = calculateresponsedistance(nKRange(indexM), nP, p, c, t);
    valStereotypy = fnStereotypy(D2, D1);
    save('data\theoretical_stereotypy.mat', 'D1', 'D2', 'nP', 'c', 'p', 'q', 't', 'nKRange', 'valStereotypy');
    disp(nKRange(indexM))
end
end % computetheoreticaltotalkcresponsestereotypy

function [D1, D2] = calculateresponsedistance(nK, nP, p, c, t)
%%CALCULATERESPONSEDISTANCE Calculates the distance between theoretical
%neuron responses for specified values

% calculate constants
v = 0:nK;
[r1, r2] = meshgrid(v, v);
inputDiff = ((r1(:) - r2(:)) .^ 2).';
v = 0:nP;
[R1, V] = meshgrid(r1(:), v);
[R2, ~] = meshgrid(r2(:), v);
R1pdf = zeros(size(R1));
R2pdf = zeros(size(R2));
tempccdf = 1 - binocdf(t - 1, v, c);
ccdf = tempccdf(V+1);
for index = 1:length(v)
    R1pdf(index, :) = binopdf(R1(index, :), nK, ccdf(index, :));
    R2pdf(index, :) = binopdf(R2(index, :), nK, ccdf(index, :));
end
tempVpdf = binopdf(v, nP, p);
Vpdf = tempVpdf(V+1);
%D1
dPdf = sum(Vpdf .* R1pdf .* R2pdf, 1);
D1 = sum(inputDiff .* dPdf, 2);
% D2
dPdf = (sum(Vpdf .* R1pdf, 1) .* sum(Vpdf .* R2pdf, 1));
D2 = sum(inputDiff .* dPdf, 2);
end % calculateresponsedistance