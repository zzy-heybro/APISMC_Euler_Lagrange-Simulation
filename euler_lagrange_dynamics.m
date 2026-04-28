function dx = euler_lagrange_dynamics(t, x,param)
global logs
    theta  = x(1:2);    dtheta = x(3:4);    xi_hat = x(5:7);  
    theta_d = [1.5 - 1.3*exp(-t) + 0.5*exp(-4*t); 1.3 + exp(-t) - 0.5*exp(-4*t)];
    dtheta_d = [1.3*exp(-t) - 2*exp(-4*t); -exp(-t) + 2*exp(-4*t)];
    ddtheta_d = [-1.3*exp(-t) + 8*exp(-4*t); exp(-t) - 8*exp(-4*t)]; e1 = theta - theta_d; de1 = dtheta - dtheta_d;
    m1 = param.m1; m2 = param.m2; l1 = param.l1; r1 = param.r1; r2 = param.r2;J1 = param.J1; J2 = param.J2; g = param.g; 
    M11 = m1*r1^2 + m2*(l1^2 + r2^2 + 2*l1*r2*cos(theta(2))) + J1 + J2; M12 = m2*r2^2 + m2*l1*r2*cos(theta(2)) + J2;  M21 = M12;  M22 = m2*r2^2 + J2;  
    C11 = -m2*l1*r2*sin(theta(2))*dtheta(2);  C12 = -m2*l1*r2*sin(theta(2))*(dtheta(1)+dtheta(2)); C21 = m2*l1*r2*sin(theta(2))*dtheta(1);   C22 = 0;    
    G1 = (m1*r1 + m2*l1)*g*cos(theta(1)) + m2*r2*g*cos(theta(1)+theta(2));    G2 = m2*r2*g*cos(theta(1)+theta(2)); 
    M = [M11 M12; M21 M22]; C = [C11 C12; C21 C22]; G = [G1; G2];  mi1=0.5;  mi2=0.3;
    Mi11 = mi1*r1^2 + mi2*(l1^2 + r2^2 + 2*l1*r2*cos(theta(2))) + J1 + J2;    Mi12 = mi2*r2^2 + mi2*l1*r2*cos(theta(2)) + J2;
    Mi21 = Mi12; Mi22 = mi2*r2^2 + J2; Mi = [Mi11 Mi12; Mi21 Mi22];    Ci11 = -mi2*l1*r2*sin(theta(2))*dtheta(2);
    Ci12 = -mi2*l1*r2*sin(theta(2))*(dtheta(1)+dtheta(2));    Ci21 = mi2*l1*r2*sin(theta(2))*dtheta(1);
    Ci22 = 0; Ci = [Ci11 Ci12; Ci21 Ci22];    Gi1 = (mi1*r1 + mi2*l1)*g*cos(theta(1)) + mi2*r2*g*cos(theta(1)+theta(2));
    Gi2 = mi2*r2*g*cos(theta(1)+theta(2));    Gi = [Gi1; Gi2]; tau_d = [3*sin(t) ; 6*cos(t) ];
    mu1 = param.mu1; mu2 = param.mu2; mu11 = param.mu11; mu22 = param.mu22;
    gam_e = gamma((1-(param.mu11 + 1)/2)/((param.mu22-param.mu11)/2)) * gamma((param.mu22 -1)/2)/((param.mu22-param.mu11)/2) * ((2^(param.mu11 + 1)/2)/((2^(param.mu22 + 1)/2)) * param.mu2 * (2^((1 - param.mu22)/2)))^ (1-(param.mu11 + 1)/2)  / 2^(param.mu11 + 1)/2 * param.mu1 * ((param.mu22-param.mu11)/2);
    gam_s = gamma((1-(param.alpha1 + 1)/2)/((param.alpha2-param.alpha1)/2)) * gamma((param.alpha2 -1)/2)/((param.alpha2-param.alpha1)/2) * ((2^(param.alpha1 + 1)/2)/((2^(param.alpha2 + 1)/2)) * param.sigma2 * (2^((1 - param.alpha2)/2)))^ (1-(param.alpha1 + 1)/2) / 2^(param.alpha1 + 1)/2 * param.sigma1 * ((param.alpha2-param.alpha1)/2); T_s = .15; T_e = .15;
    e2 = de1 + gam_e/T_e * ( mu1 * abs(e1).^mu11 .* sign(e1) + mu2 * abs(e1).^mu22 .* sign(e1));
    s = e2 + gam_s/T_s * ( param.sigma1 * integral_power(e2, param.alpha1, t) + param.sigma2 * integral_power(e2, param.alpha2, t))+  param.sigma3 * frac_diff(e2, param.beta, t).* sign(e2);
    Nval = -M \ (C * dtheta + G) - ddtheta_d;
    dxi1_dt = param.eta1 * norm(s);  dxi2_dt = param.eta2 * norm(s) * norm(theta);  dxi3_dt = param.eta3 * norm(s) * norm(dtheta)^2;
    tau = -M * (Nval + gam_e/T_e * ( param.sigma1 * abs(e2).^param.alpha1 .* sign(e2) + param.sigma2 * abs(e2).^param.alpha2 .* sign(e2)) + param.sigma3 * frac_diff(e2, param.beta, t) +  mu1*mu11*diag(abs(e1).^(mu11-1)) * de1 + mu2*mu22*diag(abs(e1).^(mu22-1)) * de1) - M * ((xi_hat(1) + xi_hat(2)*norm(theta) + xi_hat(3)*norm(dtheta)^2)*sign(s) + 1* ( param.sigma4 * abs(s).^param.gamma1 .* sign(s) +  param.sigma5 * abs(s).^param.gamma2 .* sign(s)));
    ddtheta = (M+Mi) \ (tau+tau_d- (C+Ci) * dtheta - G - Gi) ;
    dx = [dtheta; ddtheta; dxi1_dt; dxi2_dt; dxi3_dt];
    logs.t(end+1) = t; logs.e1(:,end+1) = e1; logs.s(:,end+1) = s; logs.xi_hat(:,end+1) = xi_hat;  
end
function D = frac_diff(signal, beta, t)
    persistent h N buff hist time_prev signal_dim
    if isempty(N)
        h = 0.1;
        N = 200;
        signal_dim = length(signal);
        buff = zeros(N, signal_dim);     
        hist = zeros(N, 1);              
        for k = 0:N-1
            hist(k+1) = (-1)^k * gamma(beta + 1) / (gamma(k + 1) * gamma(beta - k + 1));
        end
        time_prev = 0;
    end
    if size(signal, 2) > 1
        signal = signal(:);
    end
    if t - time_prev > 1e-9
        new_row = (abs(signal) .* sign(signal))';  
        buff = [new_row; buff(1:end-1, :)];        
        time_prev = t;
    end
    D = zeros(signal_dim, 1);
    for i = 1:signal_dim
        if length(hist) ~= size(buff(:, i), 1)
            error('hist and buff dimension mismatch.');
        end
        D(i) = sum(hist .* buff(:, i)) / h^beta;
    end
end
function I = integral_power(signal, alpha, t)
    persistent memory time_prev integral_memory
    if isempty(memory)
        memory = zeros(size(signal));
        integral_memory = zeros(size(signal));
        time_prev = 0;
    end
    dt = t - time_prev;
    if dt > 1e-6
        memory = abs(signal).^alpha .* sign(signal);
        integral_memory = integral_memory + memory * dt;
        time_prev = t;
    end
    I = integral_memory;
end