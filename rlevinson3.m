function [ ki ] = rlevinson3( Am, M )
% RLEVINSON3 Find cholesky factor, prediction error variances,
% and lattice coefficients using the reverse Levinson procedure.
%
% Inputs:
%   Am                Autocorrelation Sequence
%
% Outputs:
%   k_i               Lattice coefficients

% Curate inputs
Am = Am(:);
if ~exist('M', 'var')
    M = length(Am);
end

% Initialize outputs
ki = zeros(M,1);

% Perform recursion
for m = 1:M
    ki(m) = Am(1);
    Am = (Am + ki(m)*flipud(Am)) ./ (1-abs(ki(m))^2);
end

end

