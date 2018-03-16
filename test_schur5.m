clear all

% Demo from class
%
% sig = 3;
% alpha = 2;
% cm = @(k)(sig^2 * alpha^abs(k));
% 
% cmi = [];
% N = 20;
% for k = 0:N-1
%     cmi = [cmi; cm(k)];
% end
% 
% [ Fch, sig_m, k ] = schur5(cmi)

% Problem 5
cm = 0.5*ones(20,1); cm(1) = 1;

[Fch, sig_m, k] = schur5(cm)