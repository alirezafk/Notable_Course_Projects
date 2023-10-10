function [suc_re_rate, num_recover, Dcorrect, Scorrect] = Evaluation(D, Dhat, Shat, thre)
    [~, N] = size(D);
    recover_index = (-1)*ones(1,N);
    correct_perm = (-1)*ones(1,N);
    dir_uncertainty = ones(1, N);
    
    for i = 1:N
        valid_perm_ind = setdiff(1:N, correct_perm);
        valid_rec_ind = setdiff(1:N, recover_index);
        
        ro_perm = D(:, i)'*Dhat(:, valid_perm_ind);
        [~, idx_perm] = max(abs(ro_perm));
        correct_perm(i) = valid_perm_ind(idx_perm);
        
        if ro_perm(idx_perm) < 0
            dir_uncertainty(i) = -1;
        end
        
        ro_rec = Dhat(:, i)'*D(:, valid_rec_ind);
        [val, idx_rec] = max(abs(ro_rec));
        if val >= thre
            recover_index(i) = valid_rec_ind(idx_rec);
        end
    end
    
    num_recover = sum(recover_index > 0);
    suc_re_rate = (num_recover/N)*100;
    Dcorrect = Dhat*diag(dir_uncertainty);
    Dcorrect = Dcorrect(:, correct_perm);
    Scorrect = diag(dir_uncertainty)*Shat;
    Scorrect = Scorrect(correct_perm, :);
end