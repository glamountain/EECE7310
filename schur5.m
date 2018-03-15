function [ Fch, sig_m, k_i ] = schur5( cm )
% SCHUR5 Find cholesky factor, prediction error variances,
% and lattice coefficients using the Schur procedure.
%
% Inputs:
%   cm                Autocorrelation Sequence
%
% Outputs:
%   Fch               Cholesky factor
%   sig_m             Prediction error variances
%   k_i               Lattice coefficients

% Curate input
cm = cm(:);
cmm = cm ./ sqrt(cm(1));
G = [cmm, [0; cmm(2:end)]];

M = length(cmm);

% Initialize Outputs
Fch = zeros(M);
sig = zeros(M,1);
k_i = zeros(M,1);

for s = 1:M
    % Collect
    sig(s) = G(s,1);
    Fch(:,s) = G(:,1) ./ sig(s);
    
    % Perform shift and compute k
    G_shift = [ [ 0; G(:,1)], [ G(:,2); 0] ];
    k_i(s) = G_shift(s+1,2) / G_shift(s+1,1);
    
    % Apply transformation
    hM = [1 -k_i(s); -conj(k_i(s)), 1] * (1 - abs(k_i(s))^2)^(-0.5);
    G = G_shift * hM;
    
    G = G(1:length(cmm),:);
end

sig_m = diag(sig);
    
end

