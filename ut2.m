clear;clc
close all

%% Problem 2

N = 1e4; % Number of values of Kappa

Mx = 1;
g = @(x)(exp(x));

for a = [0.5,2]
    
    mux = 0;
    sigx2 = (a^2)/3;
    
    muy = (exp(a) - exp(-a))/(2*a);
    sigy2 = (exp(2*a) - exp(-2*a))/(4*a) - muy^2;
    
    K = zeros(1,N);
    muy_hat = zeros(1,N);
    sigy2_hat = zeros(1,N);
    
    for i = 1:N
        
        K(i) = (-Mx + i/(N/10));
        
        x1_til = sqrt((Mx + K(i))*sigx2);
        x = [mux - x1_til, mux, mux + x1_til];
        
        y = g(x);
        W = [1/(2*(Mx + K(i))), K(i)/(Mx + K(i)), 1/(2*(Mx + K(i)))];
        
        muy_hat(i) = sum(W.*y);
        sigy2_hat(i) = sum(W.*((y - muy_hat(i)).^2));
        
    end
    
    [~, opt_mu] = min(abs(muy_hat - muy));
    [~, opt_var] = min(abs(sigy2_hat - sigy2));
    
    
    set(groot, 'defaultTextInterpreter','latex');
    lw = 3;
    fs = 16;
    figure('rend','painters','pos',[10 10 900 600])
    
    subplot(2,1,1);
    plot(K, muy_hat,'LineWidth',lw); hold on
    plot(K, muy * ones(1,N), '--r','LineWidth',lw);
    title(sprintf('$\\mu_{y}$ for a=%0.2f, Best estimate at $\\kappa=%0.2f$',a,K(opt_mu)),'FontSize',fs)
    legend({'$\hat{\mu_{y}}$','$\mu_{y}$'},'Location','NorthWest','FontSize',fs)
    xlabel('$\kappa$','FontSize',fs);
    
    subplot(2,1,2);
    plot(K, sigy2_hat,'LineWidth',lw); hold on
    plot(K, sigy2 * ones(1,N), '--r','LineWidth',lw);
    title(sprintf('$\\sigma^{2}_{y}$ for a=%0.2f, Best estimate at $\\kappa=%0.2f$',a,K(opt_var)),'FontSize',fs)
    legend({'$\hat{\sigma^{2}_{y}}$','$\sigma^{2}_{y}$'},'Location','NorthWest','FontSize',fs)
    xlabel('$\kappa$','FontSize',fs);
    
end