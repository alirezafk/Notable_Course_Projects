function [Dhat, Shat] = MOD(X, N, N0, maxIter, sp_alg, thre)
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
        end
        
        % Now assume S is fixed and known
        Dhat = X*pinv(Shat);
        Dhat = Dhat./vecnorm(Dhat);
        
        curObj = X-Dhat*Shat;
        objVals(numIter) = norm(curObj(:), 2)/normX;
        if numIter > 1
            objDiff = abs(objVals(numIter) - objVals(numIter-1));
        end
        numIter = numIter + 1;
    end
    
    fprintf("MOD algorithm finished in %d iterations!\n", numIter);
    figure
    plt_range = 1:(numIter-1);
    plot(plt_range, objVals(plt_range));
    grid on
    xlabel("Iteration");
    ylabel("objective");
    title("normalized objective value for MOD: $\frac{||X-\hat{D}\hat{S}||}{||X||}$" , 'interpreter', 'latex');
end
        