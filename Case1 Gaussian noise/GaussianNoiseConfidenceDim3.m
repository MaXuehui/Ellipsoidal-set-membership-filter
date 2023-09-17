clear;
clc;

%%%%%%%%  Gaussian Noise (confidence=99%  90%) %%%%%%%%%%%%%
Dim = 3;
Sigma = 0.00072*[5 0 -1; 0 5 0; -1 0 5];    %  0.99->0.000675*[5 0 -1; 0 5 0; -1 0 5]      0.90->0.00072*[5 0 -1; 0 5 0; -1 0 5]
miu = [0; 0; 0];
Q = [5 0 -1; 0 5 0; -1 0 5];

dim1 = -3.5:0.05:3.5;   %  0.99->-3.5:0.05:3.5      0.90->-3.5:0.05:3.5
dim2 = -3.5:0.05:3.5;
dim3 = -3.5:0.05:3.5;
N = length(dim1);
dim4 = 0:0.5:350;   %  0.99->0:1:350      0.90->0:0.5:350
M = length(dim4);

index = 0;
noise_index = 0;
for i = 1:N
    for j = 1:N
        for k = 1:N
            index = index+1;
            X_dim3 = [dim1(i); dim2(j); dim3(k)];
            pdf = 1/(2*pi)^(Dim/2)*det(Sigma)^(-1/2)*exp(-1/2*(X_dim3-miu)'*Sigma^(-1)*(X_dim3-miu));
            x_plot(index) = dim1(i);
            y_plot(index) = dim2(j);
            z_plot(index) = dim3(k);
            pdf_3(index) = pdf;
            for l = 1:M
                if dim4(l)<pdf
                    noise_index = noise_index+1;
                    noise_data(:,noise_index) = [dim1(i); dim2(j); dim3(k)];
                end
            end
        end
    end
end

N_noise = 500;
random_index = randi(noise_index,1,N_noise);
for i = 1:N_noise
    random_noise_data(:,i) = noise_data(:,random_index(i));
end
figure(2)
plot3(random_noise_data(1,:),random_noise_data(2,:),random_noise_data(3,:),'.','Linewidth', 2, 'MarkerSize', 8); hold on;


syms x y z
ellipsoid_sym = [(x-0); (y-0); (z-0)]'*inv(Q)*[(x-0); (y-0); (z-0)]-1;
fimplicit3(ellipsoid_sym,'EdgeColor','none','FaceAlpha',.5); hold off;
xlabel('x');
ylabel('y');
zlabel('z');





% dim4 = 0:0.0001:0.035;
% M = length(dim3);
% noise_index = 0;
% index = 0;
% for i = 1:N
%     for j = 1:N
%         index = index+1;
%         x_dim2 = [dim1(i); dim2(j)];
%         pdf = 1/(2*pi)^(2/2)*det(Sigma)^(-1/2)*exp(-1/2*(x_dim2-miu)'*Sigma^(-1)*(x_dim2-miu));
%         x_plot(index) = dim1(i);
%         y_plot(index) = dim2(j);
%         z_plot(index) = pdf;
%         for k = 1:M
%             if dim3(k)<=pdf
%                 noise_index = noise_index+1;
%                 noise_data(:,noise_index) = [dim1(i);dim2(j)];
%             end
%         end
%     end
% end
% figure(1)
% plot3(x_plot,y_plot,z_plot,'.b');
% 
% 
% N_noise = 500;
% random_index = randi(noise_index,1,N_noise);
% for i = 1:N_noise
%     random_noise_data(:,i) = noise_data(:,random_index(i));
% end
% figure(2)
% plot(random_noise_data(1,:),random_noise_data(2,:),'.','Linewidth', 2, 'MarkerSize', 8); hold on;
% 
% 
% 
% syms x y
% ellipsoid_sym = [(x-0); (y-0)]'*inv(R)*[(x-0); (y-0)]-1;
% fimplicit(ellipsoid_sym); hold off;
% xlabel('x');
% ylabel('y');


































































