clear

%% 3b
N = 100;

alpha = 0.9;
sig_0 = 1;

f = @(x,a)([a*x; a]);
h = @(x,a)(x);

xi_true = zeros(2,N);
xi_est = zeros(size(xi_true));
y = zeros(1,N);

xi_true(:,1) = [1; alpha];
for n = 1:N
    xi_true(:,n+1) = f(xi_true(1,n),alpha) + [1;0]*randn;
    y(n) = h(xi_true(1,n),alpha);
end

%% 3c
Fk = @(x,a)([a,x; 0,1]);
Hk = @(x,a)([1,0]);
R = 0;
Q = [sig_0^2,0; 0 0];

% Initial
P_upd = eye(2);
xi_upd = [0;0];

% Recursion
for n = 1:N
    
    % Predict
    xi_pred = f(xi_upd(1),xi_upd(2));
    F = Fk(xi_upd(1),xi_upd(2));
    P_pred = F*P_upd*F' + Q;

    % Update
    y_est = y(n) - h(xi_pred(1),xi_pred(2));
    H = Hk(xi_pred(1),xi_pred(2));
    S = H*P_pred*H' + R;
    K = P_pred*H'/S;
    
    xi_upd = xi_pred + K*y_est;
    P_upd = (1 - K*H)*P_pred;
    
    xi_est(:,n) = xi_upd;
end

plot(xi_true(2,:),'LineWidth',3); hold on
plot(xi_est(2,:),'--r'); hold off
xlabel('Iteration'); ylabel('Alpha');
legend('True Value', 'Kalman Estimate');