function [x_kf_pred, P_kf_pred, yinnovations, HPH] = KF_pred(y,H,A,B,Q,x_est_ant,P_est_ant)

% KF operating on a sample basis
% 
% v0: Pau Closas (06/2016)

x_kf_pred = A*x_est_ant;
P_kf_pred = Q + A*P_est_ant*A';

HPH = H*P_kf_pred*H';
yinnovations = y - H*x_kf_pred;



