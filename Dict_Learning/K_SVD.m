function [Dhat, Shat] = K_SVD(X, N, N0, maxIter, sp_alg, thre)
    [M, ~] = size(X);
    Dhat = randn(M, N);
    Dhat = Dhat./vecnorm(Dhat);
    objDiff = inf;
    objVals = zeros(1, maxIter);
    normX = norm(X(:), 2);
    numIter = 1;
    
    while(objDiff > thre && numIter < maxIter)
        % Fist we assume D is fixed and known
        switch sp_alg
            case "MP"
                Shat = MP(X, Dhat, N0);
            case "OMP"
                Shat = OMP(X, Dhat, N0);
            otherwise
                disp("Invalid Algorithm!!");
                Dhat = nan;
                Shat = nan;
        end
        
        % Now assume S is fixed and known
        R = X - Dhat*Shat;
        for i = 1:N
            Ri = R + Dhat(:, i)*Shat(i, :);
            if sum(Shat(i,:) ~= 0) > 0
                nonzero_index = find(Shat(i,:));
                Ri_modify = Ri(:, nonzero_index);
                
                % perform svd
                [U, Sigma, V] = svd(Ri_modify);
                
                % Update
                Dhat(:, i) = U(:, 1);
                Shat(i, nonzero_index) = Sigma(1,1)*V(1,:);
            end
        end
        curObj = X-Dhat*Shat;
        objVals(numIter) = norm(curObj(:), 2)/normX;
        if numIter > 1
            objDiff = abs(objVals(numIter) - objVals(numIter-1));
        end
        numIter = numIter + 1;
    end
    
    fprintf("K-SVD algorithm finished in %d iterations!\n", numIter);
    figure
    plt_range = 1:(numIter-1);
    plot(plt_range, objVals(plt_range));
    grid on
    xlabel("Iteration");
    ylabel("objective");
    title("normalized objective value for K-SVD: $\frac{||X-\hat{D}\hat{S}||}{||X||}$" , 'interpreter', 'latex');
end