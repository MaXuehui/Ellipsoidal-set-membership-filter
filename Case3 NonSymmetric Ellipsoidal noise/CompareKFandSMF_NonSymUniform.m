clear;
clc;

%%%%%%%%  non-symmetric Noise  %%%%%%%%%%%%%
load NoiseNonSymUniformData3dim.mat
NoiseNonSymUniformData3dim = random_ellips_norm3_NonSym;
load NoiseNonSymUniformData2dim.mat
NoiseNonSymUniformData2dim = random_ellips_norm2_NonSym;

N = length(NoiseNonSymUniformData2dim);

F = [0 1 0; 0 0 1; 0.2 -0.9 1.3];
Q = [5 0 -1; 0 5 0; -1 0 5];
H = [1.2 1.5 -0.9; -1 0.8 1.1];
R = [5 -2; -2 5];

x(:,1) = [0; 0; 0];
y(:,1) = H*x(:,1)+NoiseNonSymUniformData2dim(:,1);
y_real(:,1) = H*x(:,1);

for i = 1:N-1
    x(:,i+1) = F*x(:,i)+NoiseNonSymUniformData3dim(:,i);
    y(:,i+1) = H*x(:,i+1)+NoiseNonSymUniformData2dim(:,i);
    y_real(:,i+1) = H*x(:,i+1);
end



%%%%%%%%%  KalmanFilter (confidence=99%)  %%%%%%%%%

Sigma_3Dim = 0.000675*[5 0 -1; 0 5 0; -1 0 5];
Sigma_2Dim = 0.01*[5 -2; -2 5];

x_hat(:,1) = [0; 0; 0];
y_hat(:,1) = H*x_hat(:,1);
P(:,:,1) = [10 0 0; 0 10 0; 0 0 10];
for i = 1:N-1
    x_hat(:,i+1) = F*x_hat(:,i);
    P(:,:,i+1) = F*P(:,:,i)*F'+Sigma_3Dim;
    K(:,:,i+1) = P(:,:,i+1)*H'*inv(H*P(:,:,i+1)*H'+Sigma_2Dim);
    x_hat(:,i+1) = x_hat(:,i+1)+K(:,:,i+1)*(y(i+1)-H*x_hat(:,i+1));
    P(:,:,i+1) = (eye(3,3)-K(:,:,i+1)*H)*P(:,:,i+1);
    y_hat(:,i+1) = H*x_hat(:,i+1);
end

figure(1)
plot(y(1,:)); hold on;
plot(y_real(1,:)); hold on;
plot(y_hat(1,:)); hold off;
legend('obsevation','real','estimate');
figure(2)
plot(y(2,:)); hold on;
plot(y_real(2,:)); hold on;
plot(y_hat(2,:)); hold off;
legend('obsevation','real','estimate');



%%%%%%%%%  KalmanFilter (confidence=90%)  %%%%%%%%%

Sigma_3Dim = 0.00072*[5 0 -1; 0 5 0; -1 0 5];
Sigma_2Dim = 0.2*[5 -2; -2 5];

x_hat(:,1) = [0; 0; 0];
y_hat(:,1) = H*x_hat(:,1);
P(:,:,1) = [10 0 0; 0 10 0; 0 0 10];
for i = 1:N-1
    x_hat(:,i+1) = F*x_hat(:,i);
    P(:,:,i+1) = F*P(:,:,i)*F'+Sigma_3Dim;
    K(:,:,i+1) = P(:,:,i+1)*H'*inv(H*P(:,:,i+1)*H'+Sigma_2Dim);
    x_hat(:,i+1) = x_hat(:,i+1)+K(:,:,i+1)*(y(i+1)-H*x_hat(:,i+1));
    P(:,:,i+1) = (eye(3,3)-K(:,:,i+1)*H)*P(:,:,i+1);
    y_hat(:,i+1) = H*x_hat(:,i+1);
end

figure(3)
plot(y(1,:)); hold on;
plot(y_real(1,:)); hold on;
plot(y_hat(1,:)); hold off;
legend('obsevation','real','estimate');
figure(4)
plot(y(2,:)); hold on;
plot(y_real(2,:)); hold on;
plot(y_hat(2,:)); hold off;
legend('obsevation','real','estimate');


%%%%%%%%%  set membership filter  %%%%%%%%%
P(:,:,1) = [10 0 0; 0 10 0; 0 0 10];
x_hat(:,1) = [0; 0; 0];
y_hat(:,1) = H*x_hat(:,1);
syms q_solve p_solve
for i = 1:N-1

    %%%% state update %%%%%
    x_hat(:,i+1) = F*x_hat(:,i);

    %%%%%% minimal volum p  %%%%%%%%%%%
    n = 3;
    Q1 = F*P(:,:,i)*F';
    Q2 = Q;
    eigenvalue = eig(Q1*inv(Q2));
    eqn = 1/(eigenvalue(1)+p_solve)+1/(eigenvalue(2)+p_solve)+1/(eigenvalue(3)+p_solve) == n/p_solve/(p_solve+1);
    p_ = double(solve(eqn,p_solve));
    index = find(p_>0,1);
    p = p_(index);    

    %%%%%% matrix updata %%%%%%%
    P(:,:,i+1) = (1/p+1)*F*P(:,:,i)*F'+(p+1)*Q;

    %%%%%% minimal volum p  %%%%%%%%%%%
    n = 3;
    e(:,i+1) = y(:,i+1)-H*x_hat(:,i+1);
    criteria = n*(1-e(:,i+1)'*inv(R)*e(:,i+1))-trace(P(:,:,i+1)*H'*inv(R)*H);
    if criteria > 0
        q = 0;
    else
%         rank(P(:,:,i+1)*H'*inv(R)*H)
        eigenvalue_q = eig(P(:,:,i+1)*H'*inv(R)*H);
        beta_d = 1-e(:,i+1)'*inv(q_solve^(-1)*R+H*P(:,:,i+1)*H')*q_solve^(-2)*R*inv(q_solve^(-1)*R+H*P(:,:,i+1)*H')*e(:,i+1);
        beta = 1+q_solve-e(:,i+1)'*inv(1/q_solve*R+H*P(:,:,i+1)*H')*e(:,i+1);
        eqn = eigenvalue_q(1)/(1+q_solve*eigenvalue_q(1))+eigenvalue_q(2)/(1+q_solve*eigenvalue_q(2)) == n*beta_d/beta;
%         eqn = 1/(eigenvalue(1)+q_solve)+1/(eigenvalue(2)+q_solve) == n*beta_d/beta;
        q_ = double(solve(eqn,q_solve));
        index = find(q_>0,1);
        q = q_(index);
    end

    %%%%%% obsevation update %%%%
    L(:,:,i+1) = P(:,:,i+1)*H'*inv(H*P(:,:,i+1)*H'+1/q*R);
%     e(:,i+1) = z_real(:,i+1)-H*x_hat(:,i+1);
    x_hat(:,i+1) = x_hat(:,i+1)+L(:,:,i+1)*e(:,i+1);
    beta = 1+q-e(:,i+1)'*inv(1/q*R+H*P(:,:,i+1)*H')*e(:,i+1);
    P(:,:,i+1) = beta*((eye(3,3)-L(:,:,i+1)*H)*P(:,:,i+1)*(eye(3,3)-L(:,:,i+1)*H)'+1/q*L(:,:,i+1)*R*L(:,:,i+1)');

    y_hat(:,i+1) = H*x_hat(:,i+1);

end

figure(5)
plot(y(1,:)); hold on;
plot(y_real(1,:)); hold on;
plot(y_hat(1,:)); hold off;
legend('obsevation','real','estimate');
figure(6)
plot(y(2,:)); hold on;
plot(y_real(2,:)); hold on;
plot(y_hat(2,:)); hold off;
legend('obsevation','real','estimate');
























