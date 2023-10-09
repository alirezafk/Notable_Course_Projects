function data = generate_data(modulation, N)
    switch modulation
        case "BPSK"
            data = 2*randi(2, [1,N])-3;
        case "QPSK"
            coeff1 = 2*randi(2, [1,N])-3;
            coeff2 = 2*randi(2, [1,N])-3;
            data = coeff1*(1/sqrt(2)) + coeff2*(1i/sqrt(2));
        otherwise
            disp("invalid modulation!");
            return
    end
