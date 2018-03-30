N = 100;

theta = 5;
sig_0 = 1;

f = @(x)(x);
h = @(x)(sqrt(x));

y = sqrt(theta) + sig_0*randn(1,N);
x = zeros(size(y));

Fk = 1;
Hk = @(x) 0.5 / sqrt(x);
R = sig_0^2;
Q = 0;

P_upd = 1;
x_upd = theta;
for n = 1:N
    
    % Predict
    x_pred = f(x_upd);
    P_pred = Fk*P_upd*Fk' + Q;

    % Update
    y_est = y(n) - h(x_pred);
    H = Hk(x_pred);
    S = H*P_pred*H' + R;
    K = P_pred*H'/S;
    
    x_upd = x_pred + K*y_est;
    P_upd = (1 - K*H)*P_pred;
    
    x(n) = x_upd;
end

plot(y); hold on
plot(h(x), 'r'); hold off