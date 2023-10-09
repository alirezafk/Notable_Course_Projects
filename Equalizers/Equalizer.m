function d = Equalizer(Type, ntap, x, N0, Es, ks)
    d = zeros(1,ntap);
    mid = (ntap+1)/2;
    k = (ntap-1)/2;
    L = (length(x)-1)/2;
    x1 = zeros(1, ntap);
    x1((mid - L):(mid + L)) = x;
    switch Type
        case "ZF"
            q_sym = zeros(1, ntap);
            q_sym(mid) = 1;
            d_sym = sym('d', [1,ntap]);
            eqs = sym(zeros(1, ntap));
            
            for n = -k:k
                for i = -k:k
                    if (n-i) >= -L && (n-i) <= L
                        eqs(n+k+1) = eqs(n+k+1) + d_sym(mid+i)*x1(mid+n-i);
                    end
                end
            end
            
            q_sym = solve(eqs == q_sym, d_sym);
            qFieldNames = fieldnames(q_sym);
            for i = 1:ntap
                d(i) = getfield(q_sym,qFieldNames{i});
            end
        case "MMSE"
            Ry = Es*conv(x, x);
            mid = 2*L+1;
            Ry((mid-L):(mid+L)) = Ry((mid-L):(mid+L)) + 2*N0*x;
            RIy = Es*x;
            if k > L
                R = zeros(1, 4*k+1);
                R((2*(k-L)+1):(2*(k+L)+1)) = Ry;
                Ry = R;
                
                R = zeros(1, 2*k+1);
                R((k-L+1):(k+L+1)) = RIy;
                RIy = R;
            end
            vec = 0:(2*k);
            index = vec.' - vec;
            d = Ry(index+2*k+1) \ RIy((-k:k).'+k+1).';
            d = d.';
        case "DFE"
            k1 = ks(1);
            k2 = ks(2);
            Ry = Es*conv(x, x);
            mid = 2*L+1;
            Ry((mid-L):(mid+L)) = Ry((mid-L):(mid+L)) + 2*N0*x;
            RIy = Es*x;
            if k1 > 2*L
                R = zeros(1, 2*k1+1);
                R((k1+1-2*L):(k1+1+2*L)) = Ry;
                Ry = R;
            end
            if (k1+k2) > L
                R = zeros(1, 2*(k1+k2)+1);
                R((k1+k2+1-L):(k1+k2+1+L)) = RIy;
                RIy = R;
            end
            
            % Finding coefficients
            index = (0:k1).' - (0:k1);
            mid1 = (length(Ry)+1)/2;
            R11 = Ry(index + mid1);
            index = ((-k1-1):-1).' - (0:(k2-1));
            mid2 = (length(RIy)+1)/2;
            R12 = RIy(index + mid2);
            R21 = R12';
            R22 = Es*eye(k2);
            
            % matrix A
            A = [R11, R12; R21, R22];
            b1 = RIy((-k1:0).'+mid2).';
            b2 = zeros(k2, 1);
            b = [b1; b2];
            
            % solve the system of equations Ad = b
            d = A\b;
            d = d.';
    end