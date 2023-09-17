clear
clc

Q = [5 0 -1; 0 5 0; -1 0 5];
R = [5 -2; -2 5];

syms x y z

figure(3)

ellipsoid_sym = [(x-0); (y-0)]'*inv(R)*[(x-0); (y-0)]-1;
fimplicit(ellipsoid_sym); hold on;
xlabel('x');
ylabel('y');

index = 0;
index_9 = 0;
index_1 = 0;
for i = -3:0.01:3
    for j = -3:0.01:3
        if [i; j]'*inv(R)*[i; j]<=1
            if (i>0)&&(j>0)
                index_9 = index_9+1;
                Noise_ellips_norm2_9(:,index_9) = [i; j];
            end
            if (i<0)&&(j<0)
                index_1 = index_1+1;
                Noise_ellips_norm2_1(:,index_1) = [i; j];
            end
            index = index+1;
            Noise_ellips_norm2(:,index) = [i; j];            
        end
    end
end

N_1 = 2635;
random_index_1 = randi(index_1,1,N_1);
for i = 1:N_1
    Noise_ellips_norm2_NonSym(:,i) = Noise_ellips_norm2_1(:,random_index_1(i));
end
N_9 = 23722;
for i = 1:N_9
    Noise_ellips_norm2_NonSym(:,i+N_1) = Noise_ellips_norm2_9(:,i);
end
NonSym_index = N_1+N_9;

Noise_index = 500;
random_index = randi(NonSym_index,1,Noise_index);
for i = 1:Noise_index
    random_ellips_norm2_NonSym(:,i) = Noise_ellips_norm2_NonSym(:,random_index(i));
end

plot(random_ellips_norm2_NonSym(1,:),random_ellips_norm2_NonSym(2,:),'.'); hold off;


figure(4)

ellipsoid_sym = [(x-0); (y-0); (z-0)]'*inv(Q)*[(x-0); (y-0); (z-0)]-1;
fimplicit3(ellipsoid_sym,'EdgeColor','none','FaceAlpha',.5); hold on;
xlabel('x');
ylabel('y');
zlabel('z');

index = 0;
index_dim3_9 = 0;
index_dim3_1 = 0;
for i = -3:0.02:3
    for j = -3:0.02:3
        for k = -3:0.02:3
            if [i; j; k]'*inv(Q)*[i; j; k]<=1
                if (i>0)&&(j>0)&&(k>0)
                    index_dim3_9 = index_dim3_9+1;
                    Noise_ellips_norm3_9(:,index_dim3_9) = [i; j; k];
                end
                if (i<0)&&(j<0)&&(k<0)
                    index_dim3_1 = index_dim3_1+1;
                    Noise_ellips_norm3_1(:,index_dim3_1) = [i; j; k];
                end
                index = index+1;
                Noise_ellips_norm3(:,index) = [i; j; k];                
            end
        end
    end
end
% plot3(Noise_ellips_norm(1,:),Noise_ellips_norm(2,:),Noise_ellips_norm(3,:),'.');

N_1 = 61122;
random_index_1 = randi(index_dim3_1,1,N_1);
for i = 1:N_1
    Noise_ellips_norm3_NonSym(:,i) = Noise_ellips_norm3_1(:,random_index_1(i));
end
N_9 = 550098;
for i = 1:N_9
    Noise_ellips_norm3_NonSym(:,i+N_1) = Noise_ellips_norm3_9(:,i);
end
NonSym_index = N_1+N_9;

Noise_index = 500;
random_index = randi(NonSym_index,1,Noise_index);
for i = 1:Noise_index
    random_ellips_norm3_NonSym(:,i) = Noise_ellips_norm3_NonSym(:,random_index(i));
end

plot3(random_ellips_norm3_NonSym(1,:),random_ellips_norm3_NonSym(2,:),random_ellips_norm3_NonSym(3,:),'.'); hold off;





















