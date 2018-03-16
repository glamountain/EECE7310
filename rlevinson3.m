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
ki = zeros(M-1,1);

% Perform recursion
for m = M:-1:2
    ki(m-1) = -Am(m);
    Am = (Am + ki(m-1)*flipud(Am)) ./ (1-abs(ki(m-1))^2);
    Am(end) = [];
end

end

