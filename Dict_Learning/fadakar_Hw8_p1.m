clc
clear
close all
% defining parameters
M = 3;
N = 6;
T = 1e3;
mu = 1;   % Mutual coherence

while(mu > 0.9)
    % Generating Dictionary matrix
    D = randn(M, N);

    % Normalizing columns of D
    column_norms = vecnorm(D);
    D = D * diag(1./column_norms);
    
    % obtaining mutual coherence
    A = D'*D;
    idx = eye(N);
    A = A .* (1-idx);
    A = A(~idx);
    mu = max(A);
end

%% Generating source matrix S
random_matrix = rand(N, T)*10-5;
cloc = randi(N,[T,1]);
S = zeros(T, N);
S(sub2ind([T,N],(1:(T))',cloc)) = 1;
S = permute(reshape(S,1,T,N),[2 3 1])';
S = random_matrix .* S;

%% Generating Noise matrix
Noise_std = 1e-1;
Noise = Noise_std * randn(M, T);
%% Observation Matrix
X = D*S + Noise;
%Noise = 0*Noise;
%% part a
figure
scatter3(X(1,:), X(2,:), X(3,:));
title("scatterplot of observation points");
%% MOD method
sp_alg_mod = "MP";
N0 = 1;
maxIter = 2e2;
thre = 1e-11;
cor_thre = 0.99;
suc_re_rate_mod = 0;

while(suc_re_rate_mod < 100)
    [Dhat_mod, Shat_mod] = MOD(X, N, N0, maxIter, sp_alg_mod, thre);
    [suc_re_rate_mod, num_recover_mod, Dcorrect_mod, Scorrect_mod] = Evaluation(D, Dhat_mod, Shat_mod, cor_thre);
end
fprintf("MOD method:\nnumber of recovered atoms = %d of %d\nSuccessful recovery rate = %.2f\n\n", num_recover_mod, N, suc_re_rate_mod);
%% K-SVD method
sp_alg_ksvd = "MP";
suc_re_rate_ksvd = 0;
while(suc_re_rate_ksvd < 100)
    [Dhat_ksvd, Shat_ksvd] = K_SVD(X, N, N0, maxIter, sp_alg_ksvd, thre);
    [suc_re_rate_ksvd, num_recover_ksvd, Dcorrect_ksvd, Scorrect_ksvd] = Evaluation(D, Dhat_ksvd, Shat_ksvd, cor_thre);
end
fprintf("K-SVD method:\nnumber of recovered atoms = %d of %d\nSuccessful recovery rate = %.2f\n", num_recover_ksvd, N, suc_re_rate_ksvd);
%% part d:
Sdiff_mod = Scorrect_mod - S;
E_MOD = norm(Sdiff_mod(:), 2)^2/(norm(S(:), 2)^2);

Sdiff_ksvd = Scorrect_ksvd - S;
E_KSVD = norm(Sdiff_ksvd(:), 2)^2/(norm(S(:), 2)^2);

disp("E_MOD:")
disp(E_MOD)
disp("E_KSVD:")
disp(E_KSVD)