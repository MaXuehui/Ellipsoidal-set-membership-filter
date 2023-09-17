data_90 = 0;
data_10 = 0;
R = [5 -2; -2 5];
for i = 1:noise_index
    ellipsoid_sym = [noise_data(:,i)]'*inv(R)*[noise_data(:,i)]-1;
    if ellipsoid_sym<0
        data_90 = data_90+1;
    end
end

confidence = data_90/noise_index;