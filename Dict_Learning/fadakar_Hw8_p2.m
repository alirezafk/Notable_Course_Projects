clc
clear
close all
load hw8.mat
[M, N] = size(D);
[~, T] = size(S);
actual_noise = X - D*S;
actual_obj = norm(actual_noise, 2)/norm(X(:), 2);
N0 = 3;
%% MOD method
sp_alg_mod = "OMP";
maxIter = 2e2;
thre = 1e-11;
cor_thre = 0.99;
[Dhat_mod, Shat_mod] = MOD(X, N, N0, maxIter, sp_alg_mod, thre);
[suc_re_rate_mod, num_recover_mod, Dcorrect_mod, Scorrect_mod] = Evaluation(D, Dhat_mod, Shat_mod, cor_thre);
fprintf("MOD method:\nnumber of recovered atoms = %d of %d\nSuccessful recovery rate = %.2f\n\n", num_recover_mod, N, suc_re_rate_mod);
%% K-SVD method
sp_alg_ksvd = "OMP";
[Dhat_ksvd, Shat_ksvd] = K_SVD(X, N, N0, maxIter, sp_alg_ksvd, thre);
[suc_re_rate_ksvd, num_recover_ksvd, Dcorrect_ksvd, Scorrect_ksvd] = Evaluation(D, Dhat_ksvd, Shat_ksvd, cor_thre);
fprintf("K-SVD method:\nnumber of recovered atoms = %d of %d\nSuccessful recovery rate = %.2f\n", num_recover_ksvd, N, suc_re_rate_ksvd);