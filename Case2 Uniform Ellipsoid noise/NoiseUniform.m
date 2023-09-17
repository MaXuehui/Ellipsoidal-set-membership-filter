clear
clc

Q = [5 0 -1; 0 5 0; -1 0 5];
R = [5 -2; -2 5];

syms x y z

figure(1)
%%%%%%  noises are uniformly distributed in 3D ellipsoidal set  %%%%%%%%
ellipsoid_sym = [(x-0); (y-0); (z-0)]'*inv(Q)*[(x-0); (y-0); (z-0)]-1;
fimplicit3(ellipsoid_sym,'EdgeColor','none','FaceAlpha',.5); hold on;
xlabel('x');
ylabel('y');
zlabel('z');

index = 0;
for i = -3:0.02:3
    for j = -3:0.02:3
        for k = -3:0.02:3
            if [i; j; k]'*inv(Q)*[i; j; k]<=1
                index = index+1;
                Noise_ellips_norm3(:,index) = [i; j; k];                
            end
        end
    end
end
% plot3(Noise_ellips_norm(1,:),Noise_ellips_norm(2,:),Noise_ellips_norm(3,:),'.');

N = 500;
random_index = randi(index,1,N);
for i = 1:N
    random_ellips_norm3(:,i) = Noise_ellips_norm3(:,random_index(i));
end
plot3(random_ellips_norm3(1,:),random_ellips_norm3(2,:),random_ellips_norm3(3,:),'.'); hold off;


figure(2)
%%%%%%  noises are uniformly distributed in 2D ellipsoidal set  %%%%%%%%
ellipsoid_sym = [(x-0); (y-0)]'*inv(R)*[(x-0); (y-0)]-1;
fimplicit(ellipsoid_sym); hold on;
xlabel('x');
ylabel('y');

index = 0;
for i = -3:0.02:3
    for j = -3:0.02:3
        if [i; j]'*inv(R)*[i; j]<=1
            index = index+1;
            Noise_ellips_norm2(:,index) = [i; j];            
        end
    end
end
% plot(Noise_ellips_norm2(1,:),Noise_ellips_norm2(2,:),'.')
N = 500;
random_index = randi(index,1,N);
for i = 1:N
    random_ellips_norm2(:,i) = Noise_ellips_norm2(:,random_index(i));
end
plot(random_ellips_norm2(1,:),random_ellips_norm2(2,:),'.'); hold off;


















