clear
clc

Sigma = 0.2*[5 -2; -2 5];    %  0.99->0.01*[5 -2; -2 5]      0.90->0.2*[5 -2; -2 5]
miu = [0; 0];
R = [5 -2; -2 5];

dim1 = -5:0.01:5;   %  0.99->-3:0.01:3      0.90->-5:0.01:5
dim2 = -5:0.01:5;
N = length(dim1);

% index = 0;
% for i = 1:N
%     for j = 1:N
%         index = index+1;
%         x = [dim1(i); dim2(j)];
%         pdf = 1/(2*pi)^(2/2)*det(Sigma)^(-1/2)*exp(-1/2*(x-miu)'*Sigma^(-1)*(x-miu));
%         x_plot(index) = dim1(i);
%         y_plot(index) = dim2(j);
%         z_plot(index) = pdf;
%     end
% end
% 
% plot3(x_plot,y_plot,z_plot,'.b');hold on;



dim3 = 0:0.00025:0.25;   %  0.99->0:0.001:40      0.90->0:0.00025:0.25
M = length(dim3);
noise_index = 0;
index = 0;
for i = 1:N
    for j = 1:N
        index = index+1;
        x_dim2 = [dim1(i); dim2(j)];
        pdf = 1/(2*pi)^(2/2)*det(Sigma)^(-1/2)*exp(-1/2*(x_dim2-miu)'*Sigma^(-1)*(x_dim2-miu));
        x_plot(index) = dim1(i);
        y_plot(index) = dim2(j);
        z_plot(index) = pdf;
        for k = 1:M
            if dim3(k)<pdf
                noise_index = noise_index+1;
                noise_data(:,noise_index) = [dim1(i);dim2(j)];
            end
        end
    end
end
figure(1)
plot3(x_plot,y_plot,z_plot,'.b');


N_noise = 500;
random_index = randi(noise_index,1,N_noise);
for i = 1:N_noise
    random_noise_data(:,i) = noise_data(:,random_index(i));
end
figure(2)
plot(random_noise_data(1,:),random_noise_data(2,:),'.','Linewidth', 2, 'MarkerSize', 8); hold on;



syms x y
ellipsoid_sym = [(x-0); (y-0)]'*inv(R)*[(x-0); (y-0)]-1;
fimplicit(ellipsoid_sym); hold off;
xlabel('x');
ylabel('y');


% data_90 = 0;
% data_10 = 0;
% for i = 1:noise_index
%     ellipsoid_sym = [noise_data(:,i)]'*inv(Sigma)*[noise_data(:,i)]-1;
%     if ellipsoid_sym<=0
%         data_90 = data_90+1;
%     else
%         data_10 = data_10+1;
%     end
% end
% 
% confidence = data_90/(data_90+data_10);


