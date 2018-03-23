%function [x_kf_est, x] = KF_carrier_phase(Niter, N)
close all

% signal parameters
if ~exist('Niter', 'var')
    Niter = 10;     % Monte Carlo trials
end
if ~exist('N', 'var')
    N=1000;         % number of samples
end

nx = 2;
ny = 1;
x = zeros(nx,N);
y = zeros(ny,N);

R = 2;     % true covariance
R_kf_0 = R;
sigma_v = 1;          % state covariance

T = 1;
F = [0.8 0.2; 0.1 0.7];
H = [1 0];
B=[2;1];
Q = 1;
R = 2;

x(:,1) = zeros(nx,1);

% KF parameters
x0_init_est = [0;0];
Cov0_init = [3 0; 0 1];

y(:,1) = H*x(:,1) + sqrt(R)*randn;
for ii=2:N
    x(:,ii) = F*x(:,ii-1) + B*sqrt(sigma_v)*randn;
    y(:,ii) = H*x(:,ii) + sqrt(R)*randn;
end

% method parameters (Myers)
% L = 100;        % window length

% % method parameters (Bayesian Covariance Estimation - BCE)
% mu_prior_0 = zeros(ny,1);
% kappa_prior_0 = 10;
% nu_prior_0 = 10;
% Psi_prior_0 = R*(nu_prior_0 + nx + 1);  % mode around R
% %         Psi_prior_0 = R_kf_0*(nu_prior_0 + nx + 1);  % mode around R_kf_0
%
% % reference analysis, non-informative
% mu_prior_0_ref = zeros(ny,1);
% kappa_prior_0_ref = 0;        % noninf: 0
% nu_prior_0_ref = 0;               % noninf: 0 (or -1)
% Psi_prior_0_ref = 1e2*eye(ny); %R_kf_0*(nu_prior_0_ref + nx + 1); %0;%R_kf_0*(nu_prior_0 + nx + 1);  % mode around R

% Ntransient = 100;

% init
% mu_est = zeros(ny,N);
% R_est = zeros(ny,ny,N);
x_kf_pred = zeros(size(x));
P_kf_pred = zeros(nx,nx,N);
x_kf_est = zeros(size(x));
P_kf_est = zeros(nx,nx,N);
P_kf_pred_trc = zeros(1,N);
P_kf_est_trc = zeros(1,N);
HPH = zeros(ny,ny,N);
xinnovations = zeros(size(x));
yinnovations = zeros(size(y));

for mm = 1:Niter
    
    for ii = 1:N
        % init
        x_est_ant = x0_init_est;
        P_est_ant = Cov0_init;
        %         mu_prior = mu_prior_0;
        %         kappa_prior = kappa_prior_0;        % noninf: 0
        %         nu_prior = nu_prior_0;               % noninf: -1 or 0
        %         Psi_prior = Psi_prior_0;  % mode around R
        %         mu_prior_ref = mu_prior_0_ref;
        %         kappa_prior_ref = kappa_prior_0_ref;        % noninf: 0
        %         nu_prior_ref = nu_prior_0_ref;               % noninf: -1 or 0
        %         Psi_prior_ref = Psi_prior_0_ref;  % mode around R
        R_kf = R_kf_0;
        
        % Prediction Step
        [x_kf_pred(:,ii), P_kf_pred(:,:,ii), yinnovations(:,ii), HPH(:,:,ii)] = KF_pred(y(:,ii),H,F,B,Q,x_est_ant,P_est_ant);
        P_kf_pred_trc(:,ii) = trace(P_kf_pred(:,:,ii));
        
        % Update Step
        [x_kf_est(:,ii), P_kf_est(:,:,ii), xinnovations(:,ii)] = KF_upd(yinnovations(:,ii),H,R_kf,x_kf_pred(:,ii),P_kf_pred(:,:,ii));
        P_kf_est_trc(:,ii) = trace(P_kf_est(:,:,ii));
    end
    
    % for RMSE calculation
    %     R_MSE(1,:)  = R_MSE(1,:)  + (squeeze(R_est(1,1,:) - R(1,1)).^2).';
    %     R_MSE(2,:)  = R_MSE(2,:)  + (squeeze(R_est(2,2,:) - R(2,2)).^2).';
    %     R_MSE2(1,:)  = R_MSE2(1,:)  + (squeeze(R_est2(1,1,:) - R(1,1)).^2).';
    %     R_MSE2(2,:)  = R_MSE2(2,:)  + (squeeze(R_est2(2,2,:) - R(2,2)).^2).';
    %     R_MSE3(1,:)  = R_MSE3(1,:)  + (squeeze(R_est3(1,1,:) - R(1,1)).^2).';
    %     R_MSE3(2,:)  = R_MSE3(2,:)  + (squeeze(R_est3(2,2,:) - R(2,2)).^2).';
end

% R_RMSE = sqrt(R_MSE/Niter);
% R_RMSE2 = sqrt(R_MSE2/Niter);
% R_RMSE3 = sqrt(R_MSE3/Niter);

figure,
plot(x(1,:),x(2,:),'LineWidth',2), hold on;%, plot(y(1,:),'ok'),
plot(x_kf_est(1,:),x_kf_est(2,:),'r')
xlabel('time [s]'), ylabel('position [m]'), grid on
%legend('true position evolution','noisy measurements','KF (BCE)')
legend('true position evolution','KF')
% export_fig KF_est -eps -transparent

figure,
plot(P_kf_pred_trc(1,:),'LineWidth',2), hold on
plot(P_kf_est_trc(1,:),'r');
xlabel('time [s]'), ylabel('position [m]'), grid on
%legend('true position evolution','noisy measurements','KF (BCE)')
legend('Predicted Error Covariance','Updated Error Covariance')
% export_fig KF_est -eps -transparent

%end