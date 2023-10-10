Wireless Communication
Computer Assignment 3
Seyyed Mohammad Matin AleMohammad
810197457

Part1: Narrowband Channel
1
clc
clear
close all

M = 1e6;
a = 1;
x = randi([0, 1], M, 1);
x(x==0) = -1;

snr_dB = -20:1:20;
err_prob = zeros(3, length(snr_dB));
for i = 1:length(err_prob)
    h1 = sqrt(1/2)*randn(M, 1) + 1i*sqrt(1/2)*randn(M, 1);
    h2 = ones(M, 1);
    h = [h1, h2];

    snr = 10^(0.1*snr_dB(i));
    N0 = a^2/snr;
    w = sqrt(N0/2)*randn(M, 1) + 1i*sqrt(N0/2)*randn(M, 1);
    y = (h.*x + w);
    r = 2*(real(y)>0)-1;

    err_prob(1:2, i) = mean(~(r==x), 1);
    err_prob(3, i) = qfunc(sqrt(2*snr));
end

semilogy(snr_dB, err_prob);
title('Optimal Bit Error Probability for BPSK')
xlabel('SNR [dB]')
ylabel('Bit Error Probability')
legend('Fading', 'Simulation', 'Theory', 'Location', 'southwest')
disp(10*log10(qfuncinv(1e-6)^2/2))

2
clc
close all

M = 1e6;
a = sqrt(2);
Es_avg = a^2/2;
x = randi([0, 1], M, 2);
x(:, 2) = 1 - x(:, 1);
x = a*x;

snr_dB = -20:1:20;
err_prob2 = zeros(2, length(snr_dB));
for i = 1:length(err_prob2)
    snr = 10^(0.1*snr_dB(i));
    N0 = Es_avg/snr;
    w = sqrt(N0/2)*randn(M, 2) + 1i*sqrt(N0/2)*randn(M, 2);
    h = sqrt(1/2)*randn(M, 2) + 1i*sqrt(1/2)*randn(M, 2);

    y = (h.*x + w);
    [~, r] = max(abs(y).^2, [], 2);
    r = a*(r-1);
    r = (r==x(:, 1));

    err_prob2(1, i) = mean(r);
    err_prob2(2, i) = 1./(2+2*snr);
end

semilogy(snr_dB, err_prob2);
title('Optimal Bit Error Probability for BPSK')
xlabel('SNR [dB]')
ylabel('Bit Error Probability')
legend('Simulation', 'Theory', 'Location', 'southwest')
disp(10*log10(1/2/1e-6-1))

3
clc
close all

M = 1e6;
a = 1;
x = randi([0, 1], M, 1);
x(x==0) = -1;

snr_dB = -20:1:20;
err_prob3 = zeros(2, length(snr_dB));
for i = 1:length(err_prob3)
    snr = 10^(0.1*snr_dB(i));
    N0 = a^2/snr;
    w = sqrt(N0/2)*randn(M, 1) + 1i*sqrt(N0/2)*randn(M, 1);
    h = sqrt(1/2)*randn(M, 1) + 1i*sqrt(1/2)*randn(M, 1);

    y = (h.*x + w).*conj(h./abs(h));
    r = 2*(real(y)>0)-1;

    err_prob3(1, i) = mean(~(r==x), 1);
    err_prob3(2, i) = 0.5*(1-sqrt(snr/(1+snr)));
end

semilogy(snr_dB, err_prob3);
title('Optimal Bit Error Probability for BPSK')
xlabel('SNR [dB]')
ylabel('Bit Error Probability')
legend('Simulation', 'Theroty', 'Location', 'southwest')
disp(10*log10((1-2*1e-6)^2/(1-(1-2*1e-6)^2)))
semilogy(snr_dB, err_prob2(1, :), snr_dB, err_prob3(1, :));
title('Optimal Bit Error Probability for BPSK')
xlabel('SNR [dB]')
ylabel('Bit Error Probability')
legend('CSI not available at Rx but coding', 'CSI available at Rx', 'Location', 'southwest');

4
clc
close all

M = 1e6;
a = 1;
x = randi([0, 1], M, 2);
x(x==0) = -1;
x = sqrt(1/2)*(x(:, 1)+1i*x(:, 2));

snr_dB = -20:1:20;
err_prob4 = zeros(2, length(snr_dB));
for i = 1:length(err_prob4)
    snr = 10^(0.1*snr_dB(i));
    N0 = a^2/snr;
    w = sqrt(N0/2)*randn(M, 1) + 1i*sqrt(N0/2)*randn(M, 1);
    h = sqrt(1/2)*randn(M, 1) + 1i*sqrt(1/2)*randn(M, 1);

    y = (h.*x + w).*conj(h./abs(h));
    re = sqrt(1/2)*(2*(real(y)>0)-1);
    im = sqrt(1/2)*(2*(imag(y)>0)-1);

    err_prob4(1, i) = 1-(mean(re==real(x))+mean(im==imag(x)))/2;
    err_prob4(2, i) = 0.5*(1-sqrt(snr/(2+snr)));
end

semilogy(snr_dB, err_prob4);
title('Optimal Bit Error Probability for QPSK')
xlabel('SNR [dB]')
ylabel('Bit Error Probability')
legend('Simulation', 'Theory', 'Location', 'southwest')
semilogy(snr_dB, err_prob4(1,:), snr_dB, err_prob2(1,:));
title('Optimal Bit Error Probability')
xlabel('SNR [dB]')
ylabel('Bit Error Probability')
legend('QPSK', 'Orthogonal Coding', 'Location', 'southwest')

clc
close all

M = 2e6;
a = 1;
snr_dB = -10:1:10;
err_prob5_t = zeros(5, length(snr_dB));
err_prob5_s = zeros(5, length(snr_dB));
% figure;
for L = 2:2
    x = randi([0, 1], M, 1);
    x(x==0) = -1;
    x = repmat(x, 1, L);

    for i = 1:length(snr_dB)
        snr = 10^(0.1*snr_dB(i));
        N0 = a^2/snr;
        w = sqrt(N0/2)*randn(M, L) + 1i*sqrt(N0/2)*randn(M, L);
        h = sqrt(1/2)*randn(M, L) + 1i*sqrt(1/2)*randn(M, L);

        y = sum((h.*x + w).*conj(h./sqrt(sum(conj(h).*h))), 2);
        r = 2*(real(y)>0)-1;

        err_prob5_s(L, i) = mean(~(r==x(:, 1)), 1);
        coeff = 0;
        mu = sqrt(snr/(1+snr));
        for l = (1:L)-1
            coeff = coeff + nchoosek(L-1+l, l)*(1/2+mu/2)^l;
        end
        err_prob5_t(L, i) = (1/2-mu/2)^L*coeff;
        %         err_prob5_t(L, i) = nchoosek(2*L-1, L)*(1/snr/4)^L;

    end
end
semilogy(snr_dB, err_prob5_s);
title('Optimal Bit Error Probability for BPSK')
xlabel('SNR [dB]')
ylabel('Bit Error Probability')
hold on
semilogy(snr_dB, err_prob5_t, '*k');
legend('L=1', 'L=2', 'L=3', 'L=4', 'L=5', 'Theory', 'Location', 'southwest')

6
clc
close all

M = 1e6;
a = 1;
L = 2;
snr_dB = -10:1:10;
err_prob6 = zeros(2, length(snr_dB));

x = randi([0, 1], M, 1);
x(x==0) = -1;
temp1 = zeros(size(x));    temp1(2:2:end) = conj(x(1:2:end));
temp2 = zeros(size(x));    temp2(1:2:end) = -conj(x(2:2:end));
x = [x, temp1+temp2];

for i = 1:length(snr_dB)
    snr = 10^(0.1*snr_dB(i));
    N0 = a^2/snr;
    w = sqrt(N0/2)*randn(M, 1) + 1i*sqrt(N0/2)*randn(M, 1);
    h = sqrt(1/2)*randn(M, L) + 1i*sqrt(1/2)*randn(M, L);
    h(2:2:end, :) = h(1:2:end, :);

    y = sum(h.*x, 2) + w;
    y = reshape(y, 2, []);
    y(2, :) = conj(y(2, :));

    h = h(1:2:end, :).';
    y1 = sum([conj(h(1, :)); h(2, :)].*y, 1);
    y2 = sum([conj(h(2, :)); -h(1, :)].*y, 1);

    r = 2*(real([y1; y2].')>0)-1;

    err_prob6(1, i) = mean(~(r==x(1:2:end, :)), 'all');
    coeff = 0;
    mu = sqrt(snr/(1+snr));
    for l = (1:L)-1
        coeff = coeff + nchoosek(L-1+l, l)*(1/2+mu/2)^l;
    end
    err_prob6(L, i) = (1/2-mu/2)^L*coeff;

end

semilogy(snr_dB, err_prob6);
title('Optimal Bit Error Probability for Alamouti Coding')
xlabel('SNR [dB]')
ylabel('Bit Error Probability')
legend('Simulation', 'Theory', 'Location', 'southwest')

Part2: Wideband Channel
5.
clc
clear
close all

W = 20e6;
Td = 10e-6;
Tc = 5e-3;
L = fix(Td*W);
cp = L-1;
nc = 2000;
N = 1e6;
N0 = 1;

d = randi([0, 1], N, 1);
d(d==0) = -1;
d = reshape(d, nc, []);

h = sqrt(1/2)*randn(L, 1) + 1i*sqrt(1/2)*randn(L, 1);
h = [h; zeros(nc-length(h), 1)];
H = fft(h);

w = sqrt(N0/2)*randn(nc, 1) + 1i*sqrt(N0/2)*randn(nc, 1);

snr_dB = -40:1:40;
channel_capacity = zeros(size(snr_dB));
for i = 1:length(snr_dB)
    snr = 10^(0.1*snr_dB(i));
    Pmax = nc*N0*snr;

    func = @(lambda) sum(max(1/lambda-N0./abs(H).^2, 0))-Pmax;
    lambda = fzero(func, [1e-6, 1/min(N0./abs(H).^2)]); %1/max(N0./abs(H).^2)
    P = max(1/lambda-N0./abs(H).^2, 0);
    channel_capacity(i) = sum(log10(1+P.*abs(H).^2/N0));

end

plot(snr_dB, channel_capacity);
title('Channel Capacity in OFDM with Waterfilling')
ylabel('Channel Capacity [bits/OFDM Symbols]')
xlabel('SNR [dB]')

snr_dB = 0:1:60;
err_prob5_2 = zeros(size(snr_dB));
for i = 1:length(snr_dB)
    snr = 10^(0.1*snr_dB(i));
    Pmax = nc*N0*snr;

    err_prob = 0;
    for j = 1:size(d, 2)
        h = sqrt(1/2)*randn(L, 1) + 1i*sqrt(1/2)*randn(L, 1);
        h = [h; zeros(nc-length(h), 1)];
        H = fft(h);
        w = sqrt(N0/2)*randn(nc, 1) + 1i*sqrt(N0/2)*randn(nc, 1);

        func = @(lambda) sum(max(1/lambda-N0./abs(H).^2, 0))-Pmax;
        lambda = fzero(func, [1e-8, 1/min(N0./abs(H).^2)]); %1/max(N0./abs(H).^2)
        P = max(1/lambda-N0./abs(H).^2, 0);
        %         X = sqrt(P).*d(:, j);
        X = sqrt(P).*d(:, j).*conj(H./abs(H));
        x = ifft(X);
        %         x = [x(nc-cp+1:end); x];
        y = cconv(h, x, nc) + w;
        %         y = y.*conj(h./abs(h));
        Y = fft(y);

        r = 2*(real(Y)>0)-1;
        err_prob = err_prob + mean(~(r==d(:, j)));
    end

    err_prob5_2(i) = err_prob/size(d, 2);
end

semilogy(snr_dB, err_prob5_2);
title('Bit Error Probability in OFDM with Waterfilling')
xlabel('SNR [dB]')
ylabel('Bit Error Probability')

6.
clc
close all

W = 20e6;
Td = 10e-6;
Tc = 5e-3;
L = fix(Td*W);
cp = L-1;
nc = 2000;
N = 1e6;
N0 = 1;
M = 10;

d = randi([0, 1], N, 1);
d(d==0) = -1;
d = reshape(d, nc, []);

snr_dB = 0:1:20;
err_prob6_2 = zeros(1, length(snr_dB));
for i = 1:length(snr_dB)
    snr = 10^(0.1*snr_dB(i));
    Pmax = nc*N0*snr;

    err_prob = 0;
    for j = 1:size(d, 2)

        Y_MRC = zeros(nc, 1);
        for k = 1:M
            h = sqrt(1/2)*randn(L, 1) + 1i*sqrt(1/2)*randn(L, 1);
            h = [h; zeros(nc-length(h), 1)];
            H = fft(h);
            w = sqrt(N0/2)*randn(nc, 1) + 1i*sqrt(N0/2)*randn(nc, 1);

            P = Pmax/nc;
            X = sqrt(P).*d(:, j);
            x = ifft(X);
            y = cconv(h, x, nc) + w;
            Y = fft(y);
            Y = Y.*conj(H);
            Y_MRC = Y_MRC+Y;
        end

        r = 2*(real(Y_MRC)>0)-1;
        err_prob = err_prob + mean(~(r==d(:, j)));
    end
    
    err_prob6_2(i) = err_prob/size(d, 2);
end
semilogy(snr_dB, err_prob6_2);
title('Bit Error Probability in OFDM with MRC in Rx')
xlabel('SNR [dB]')
ylabel('Bit Error Probability')

7.
clc
close all

W = 20e6;
Td = 10e-6;
Tc = 5e-3;
L = fix(Td*W);
cp = L-1;
nc = 2000;
N = 1e6;
N0 = 1;

d = randi([0, 1], N, 1);
d(d==0) = -1;
d = reshape(d, nc, []);

snr_dB = 0:1:60;
err_prob7_2 = zeros(2, length(snr_dB));
for i = 1:length(snr_dB)
    snr = 10^(0.1*snr_dB(i));
    Pmax = nc*N0*snr;

    err_prob = [0; 0];
    for j = 1:size(d, 2)

        h = sqrt(1/2)*randn(L, 1) + 1i*sqrt(1/2)*randn(L, 1);
        h = [h; zeros(nc-length(h), 1)];
        H = fft(h);
        w = sqrt(N0/2)*randn(nc, 1) + 1i*sqrt(N0/2)*randn(nc, 1);

        P = Pmax/nc;
        X = sqrt(P).*d(:, j);
        x = ifft(X);
        y = cconv(h, x, nc) + w;
        Y = fft(y);
        Y_ZF = Y./H;
        Y_MMSE = Y.*conj(H./(abs(H).^2+N0*nc/Pmax));

        r_ZF = 2*(real(Y_ZF)>0)-1;
        r_MMSE = 2*(real(Y_MMSE)>0)-1;
        err_prob(1) = err_prob(1) + mean(~(r_ZF==d(:, j)));
        err_prob(2) = err_prob(2) + mean(~(r_MMSE==d(:, j)));
    end

    err_prob7_2(:, i) = err_prob/size(d, 2);
end
semilogy(snr_dB, err_prob7_2(1, :), '-*', snr_dB, err_prob7_2(2, :), '-');
title('Bit Error Probability in OFDM with Equalizer')
xlabel('SNR [dB]')
ylabel('Bit Error Probability')
legend('ZF', 'MMSE', 'Location', 'southwest');

8.
clc
close all

W = 20e6;
Td = 10e-6;
Tc = 5e-3;
L = fix(Td*W);
cp = L-1;
nc = 2000;
N = 1e6;
N0 = 1;

d = randi([0, 1], N, 1);
d(d==0) = -1;
d = reshape(d, nc, []);
% Narrowband Channel (flat fading)
snr_dB = 0:1:60;
err_prob8_2 = zeros(1, length(snr_dB));
for i = 1:length(snr_dB)
    snr = 10^(0.1*snr_dB(i));
    Pmax = nc*N0*snr;

    err_prob = 0;
    for j = 1:size(d, 2)

        h = sqrt(1/2)*randn(L, 1) + 1i*sqrt(1/2)*randn(L, 1);
%         h = [h; zeros(nc-length(h), 1)];
        H = fft(h, nc);
        w = sqrt(N0/2)*randn(nc, 1) + 1i*sqrt(N0/2)*randn(nc, 1);

        P = Pmax/nc;
        X = sqrt(P).*d(:, j);
        x = ifft(X);
        x_ceil = 0.8*max(abs(x));
        x = (min(x_ceil, abs(x))./abs(x)).*x;
        y = cconv(h, x, nc) + w;
        Y = fft(y);
        Y_MMSE = Y.*conj(H./(abs(H).^2+N0*nc/Pmax));

        r_MMSE = 2*(real(Y_MMSE)>0)-1;
        err_prob = err_prob + sum(~(r_MMSE==d(:, j)));
    end
    
    err_prob8_2(:, i) = err_prob/N;
end
semilogy(snr_dB, err_prob8_2);
title('Bit Error Probability in OFDM with MMSE and Clipping')
xlabel('SNR [dB]')
ylabel('Bit Error Probability')

semilogy(snr_dB, err_prob5_2, '-', snr_dB(1:21), err_prob6_2, '-*', snr_dB, err_prob7_2(1, :), '-o', snr_dB, err_prob7_2(2, :), '-+', snr_dB, err_prob8_2, '-d');
title('Bit Error Probability in OFDM')
xlabel('SNR [dB]')
ylabel('Bit Error Probability')
legend('Waterfilling', 'MRC(Rx)', 'ZF', 'MMSE', 'Clipping(MMSE)', 'Location', 'northeast')

% Part 3
clc
close all

M = 1e6;
a = 1;
x = randi([0, 1], M, 1);
x(x==0) = -1;

snr_dB = -20:1:20;
err_prob3 = zeros(2, length(snr_dB));
for i = 1:length(err_prob3)
    snr = 10^(0.1*snr_dB(i));
    N0 = a^2/snr;
    w = sqrt(N0/2)*randn(M, 1) + 1i*sqrt(N0/2)*randn(M, 1);
    h = sqrt(1/2)*randn(M, 1) + 1i*sqrt(1/2)*randn(M, 1);

    y = (h.*x + w).*conj(h./abs(h));
    r = 2*(real(y)>0)-1;

    err_prob3(1, i) = mean(~(r==x), 1);
    err_prob3(2, i) = 0.5*(1-sqrt(snr/(1+snr)));
end

semilogy(snr_dB, err_prob3);
title('Optimal Bit Error Probability for BPSK')
xlabel('SNR [dB]')
ylabel('Bit Error Probability')
legend('Simulation', 'Theroty', 'Location', 'southwest')

semilogy(snr_dB, err_prob2(1, :), snr_dB, err_prob3(1, :));
title('Optimal Bit Error Probability for BPSK')
xlabel('SNR [dB]')
ylabel('Bit Error Probability')
legend('CSI not available at Rx but coding', 'CSI available at Rx', 'Location', 'southwest');

% Part 4
clc
close all

M = 1e6;
a = 1;
x = randi([0, 1], M, 2);
x(x==0) = -1;
x = sqrt(1/2)*(x(:, 1)+1i*x(:, 2));

snr_dB = -20:1:20;
err_prob4 = zeros(2, length(snr_dB));
for i = 1:length(err_prob4)
    snr = 10^(0.1*snr_dB(i));
    N0 = a^2/snr;
    w = sqrt(N0/2)*randn(M, 1) + 1i*sqrt(N0/2)*randn(M, 1);
    h = sqrt(1/2)*randn(M, 1) + 1i*sqrt(1/2)*randn(M, 1);

    y = (h.*x + w).*conj(h./abs(h));
    re = sqrt(1/2)*(2*(real(y)>0)-1);
    im = sqrt(1/2)*(2*(imag(y)>0)-1);

    err_prob4(1, i) = 1-(mean(re==real(x))+mean(im==imag(x)))/2;
    err_prob4(2, i) = 0.5*(1-sqrt(snr/(2+snr)));
end

semilogy(snr_dB, err_prob4);
title('Optimal Bit Error Probability for QPSK')
xlabel('SNR [dB]')
ylabel('Bit Error Probability')
legend('Simulation', 'Theory', 'Location', 'southwest')


clc
close all

M = 2e6;
a = 1;
snr_dB = -10:1:10;
err_prob5_t = zeros(5, length(snr_dB));
err_prob5_s = zeros(5, length(snr_dB));
figure;
for L = 2:2
    x = randi([0, 1], M, 1);
    x(x==0) = -1;
    x = repmat(x, 1, L);

    for i = 1:length(snr_dB)
        snr = 10^(0.1*snr_dB(i));
        N0 = a^2/snr;
        w = sqrt(N0/2)*randn(M, L) + 1i*sqrt(N0/2)*randn(M, L);
        h = sqrt(1/2)*randn(M, L) + 1i*sqrt(1/2)*randn(M, L);

        y = sum((h.*x + w).*conj(h./sqrt(sum(conj(h).*h))), 2);
        r = 2*(real(y)>0)-1;

        err_prob5_s(L, i) = mean(~(r==x(:, 1)), 1);
        coeff = 0;
        mu = sqrt(snr/(1+snr));
        for l = (1:L)-1
            coeff = coeff + nchoosek(L-1+l, l)*(1/2+mu/2)^l;
        end
        err_prob5_t(L, i) = (1/2-mu/2)^L*coeff;
        %         err_prob5_t(L, i) = nchoosek(2*L-1, L)*(1/snr/4)^L;

    end
end
semilogy(snr_dB, err_prob5_s);
title('Optimal Bit Error Probability for BPSK')
xlabel('SNR [dB]')
ylabel('Bit Error Probability')
hold on
semilogy(snr_dB, err_prob5_t, '*k');
legend('L=1', 'L=2', 'L=3', 'L=4', 'L=5', 'Theory', 'Location', 'southwest')


% Part 6
clc
close all

M = 1e6;
a = 1;
L = 2;
snr_dB = -10:1:10;
err_prob6 = zeros(2, length(snr_dB));

x = randi([0, 1], M, 1);
x(x==0) = -1;
temp1 = zeros(size(x));    temp1(2:2:end) = conj(x(1:2:end));
temp2 = zeros(size(x));    temp2(1:2:end) = -conj(x(2:2:end));
x = [x, temp1+temp2];

for i = 1:length(snr_dB)
    snr = 10^(0.1*snr_dB(i));
    N0 = a^2/snr;
    w = sqrt(N0/2)*randn(M, 1) + 1i*sqrt(N0/2)*randn(M, 1);
    h = sqrt(1/2)*randn(M, L) + 1i*sqrt(1/2)*randn(M, L);
    h(2:2:end, :) = h(1:2:end, :);

    y = sum(h.*x, 2) + w;
    y = reshape(y, 2, []);
    y(2, :) = conj(y(2, :));

    h = h(1:2:end, :).';
    y1 = sum([conj(h(1, :)); h(2, :)].*y, 1);
    y2 = sum([conj(h(2, :)); -h(1, :)].*y, 1);

    r = 2*(real([y1; y2].')>0)-1;

    err_prob6(1, i) = mean(~(r==x(1:2:end, :)), 'all');
    coeff = 0;
    mu = sqrt(snr/(1+snr));
    for l = (1:L)-1
        coeff = coeff + nchoosek(L-1+l, l)*(1/2+mu/2)^l;
    end
    err_prob6(L, i) = (1/2-mu/2)^L*coeff;

end

semilogy(snr_dB, err_prob6);
title('Optimal Bit Error Probability for Alamouti Coding')
xlabel('SNR [dB]')
ylabel('Bit Error Probability')
legend('Simulation', 'Theory', 'Location', 'southwest')

%%%% Wide Band Channel Part 1 %%%%%%
clc
clear
close all

W = 20e6;
Td = 10e-6;
Tc = 5e-3;
L = fix(Td*W);
cp = L-1;
nc = 2000;
N = 1e6;
N0 = 1;

d = randi([0, 1], N, 1);
d(d==0) = -1;
d = reshape(d, nc, []);

h = sqrt(1/2)*randn(L, 1) + 1i*sqrt(1/2)*randn(L, 1);
h = [h; zeros(nc-length(h), 1)];
H = fft(h);

w = sqrt(N0/2)*randn(nc, 1) + 1i*sqrt(N0/2)*randn(nc, 1);

snr_dB = -40:1:40;
channel_capacity = zeros(size(snr_dB));
for i = 1:length(snr_dB)
    snr = 10^(0.1*snr_dB(i));
    Pmax = nc*N0*snr;

    func = @(lambda) sum(max(1/lambda-N0./abs(H).^2, 0))-Pmax;
    lambda = fzero(func, [1e-6, 1/min(N0./abs(H).^2)]); %1/max(N0./abs(H).^2)
    P = max(1/lambda-N0./abs(H).^2, 0);
    channel_capacity(i) = sum(log10(1+P.*abs(H).^2/N0));

end

plot(snr_dB, channel_capacity);
title('Channel Capacity in OFDM with Waterfilling')
ylabel('Channel Capacity [bits/OFDM Symbols]')
xlabel('SNR [dB]')

snr_dB = 0:1:60;
err_prob5_2 = zeros(size(snr_dB));
for i = 1:length(snr_dB)
    snr = 10^(0.1*snr_dB(i));
    Pmax = nc*N0*snr;

    err_prob = 0;
    for j = 1:size(d, 2)
        h = sqrt(1/2)*randn(L, 1) + 1i*sqrt(1/2)*randn(L, 1);
        h = [h; zeros(nc-length(h), 1)];
        H = fft(h);
        w = sqrt(N0/2)*randn(nc, 1) + 1i*sqrt(N0/2)*randn(nc, 1);

        func = @(lambda) sum(max(1/lambda-N0./abs(H).^2, 0))-Pmax;
        lambda = fzero(func, [1e-8, 1/min(N0./abs(H).^2)]); %1/max(N0./abs(H).^2)
        P = max(1/lambda-N0./abs(H).^2, 0);
        %         X = sqrt(P).*d(:, j);
        X = sqrt(P).*d(:, j).*conj(H./abs(H));
        x = ifft(X);
        %         x = [x(nc-cp+1:end); x];
        y = cconv(h, x, nc) + w;
        %         y = y.*conj(h./abs(h));
        Y = fft(y);

        r = 2*(real(Y)>0)-1;
        err_prob = err_prob + mean(~(r==d(:, j)));
    end

    err_prob5_2(i) = err_prob/size(d, 2);
end

semilogy(snr_dB, err_prob5_2);
title('Bit Error Probability in OFDM with Waterfilling')
xlabel('SNR [dB]')
ylabel('Bit Error Probability')

% Part 6
clc
close all

W = 20e6;
Td = 10e-6;
Tc = 5e-3;
L = fix(Td*W);
cp = L-1;
nc = 2000;
N = 1e6;
N0 = 1;
M = 10;

d = randi([0, 1], N, 1);
d(d==0) = -1;
d = reshape(d, nc, []);

snr_dB = 0:1:20;
err_prob6_2 = zeros(1, length(snr_dB));
for i = 1:length(snr_dB)
    snr = 10^(0.1*snr_dB(i));
    Pmax = nc*N0*snr;

    err_prob = 0;
    for j = 1:size(d, 2)

        Y_MRC = zeros(nc, 1);
        for k = 1:M
            h = sqrt(1/2)*randn(L, 1) + 1i*sqrt(1/2)*randn(L, 1);
            h = [h; zeros(nc-length(h), 1)];
            H = fft(h);
            w = sqrt(N0/2)*randn(nc, 1) + 1i*sqrt(N0/2)*randn(nc, 1);

            P = Pmax/nc;
            X = sqrt(P).*d(:, j);
            x = ifft(X);
            y = cconv(h, x, nc) + w;
            Y = fft(y);
            Y = Y.*conj(H);
            Y_MRC = Y_MRC+Y;
        end

        r = 2*(real(Y_MRC)>0)-1;
        err_prob = err_prob + mean(~(r==d(:, j)));
    end
    
    err_prob6_2(i) = err_prob/size(d, 2);
end
semilogy(snr_dB, err_prob6_2);
title('Bit Error Probability in OFDM with MRC in Rx')
xlabel('SNR [dB]')
ylabel('Bit Error Probability')

% Part 7
clc
close all

W = 20e6;
Td = 10e-6;
Tc = 5e-3;
L = fix(Td*W);
cp = L-1;
nc = 2000;
N = 1e6;
N0 = 1;

d = randi([0, 1], N, 1);
d(d==0) = -1;
d = reshape(d, nc, []);

snr_dB = 0:1:60;
err_prob7_2 = zeros(2, length(snr_dB));
for i = 1:length(snr_dB)
    snr = 10^(0.1*snr_dB(i));
    Pmax = nc*N0*snr;

    err_prob = [0; 0];
    for j = 1:size(d, 2)

        h = sqrt(1/2)*randn(L, 1) + 1i*sqrt(1/2)*randn(L, 1);
        h = [h; zeros(nc-length(h), 1)];
        H = fft(h);
        w = sqrt(N0/2)*randn(nc, 1) + 1i*sqrt(N0/2)*randn(nc, 1);

        P = Pmax/nc;
        X = sqrt(P).*d(:, j);
        x = ifft(X);
        y = cconv(h, x, nc) + w;
        Y = fft(y);
        Y_ZF = Y./H;
        Y_MMSE = Y.*conj(H./(abs(H).^2+N0*nc/Pmax));

        r_ZF = 2*(real(Y_ZF)>0)-1;
        r_MMSE = 2*(real(Y_MMSE)>0)-1;
        err_prob(1) = err_prob(1) + mean(~(r_ZF==d(:, j)));
        err_prob(2) = err_prob(2) + mean(~(r_MMSE==d(:, j)));
    end

    err_prob7_2(:, i) = err_prob/size(d, 2);
end
semilogy(snr_dB, err_prob7_2(1, :), '-*', snr_dB, err_prob7_2(2, :), '-');
title('Bit Error Probability in OFDM with Equalizer')
xlabel('SNR [dB]')
ylabel('Bit Error Probability')
legend('ZF', 'MMSE', 'Location', 'southwest');


% Part 8
clc
close all

W = 20e6;
Td = 10e-6;
Tc = 5e-3;
L = fix(Td*W);
cp = L-1;
nc = 2000;
N = 1e6;
N0 = 1;

d = randi([0, 1], N, 1);
d(d==0) = -1;
d = reshape(d, nc, []);

snr_dB = 0:1:60;
err_prob8_2 = zeros(1, length(snr_dB));
for i = 1:length(snr_dB)
    snr = 10^(0.1*snr_dB(i));
    Pmax = nc*N0*snr;

    err_prob = 0;
    for j = 1:size(d, 2)

        h = sqrt(1/2)*randn(L, 1) + 1i*sqrt(1/2)*randn(L, 1);
%         h = [h; zeros(nc-length(h), 1)];
        H = fft(h, nc);
        w = sqrt(N0/2)*randn(nc, 1) + 1i*sqrt(N0/2)*randn(nc, 1);

        P = Pmax/nc;
        X = sqrt(P).*d(:, j);
        x = ifft(X);
        x_ceil = 0.8*max(abs(x));
        x = (min(x_ceil, abs(x))./abs(x)).*x;
        y = cconv(h, x, nc) + w;
        Y = fft(y);
        Y_MMSE = Y.*conj(H./(abs(H).^2+N0*nc/Pmax));

        r_MMSE = 2*(real(Y_MMSE)>0)-1;
        err_prob = err_prob + sum(~(r_MMSE==d(:, j)));
    end
    
    err_prob8_2(:, i) = err_prob/N;
end
semilogy(snr_dB, err_prob8_2);
title('Bit Error Probability in OFDM with MMSE and Clipping')
xlabel('SNR [dB]')
ylabel('Bit Error Probability')


semilogy(snr_dB, err_prob5_2, '-', snr_dB(1:21), err_prob6_2, '-*', snr_dB, err_prob7_2(1, :), '-o', snr_dB, err_prob7_2(2, :), '-+', snr_dB, err_prob8_2, '-d');
title('Bit Error Probability in OFDM')
xlabel('SNR [dB]')
ylabel('Bit Error Probability')
legend('Waterfilling', 'MRC(Rx)', 'ZF', 'MMSE', 'Clipping(MMSE)', 'Location', 'northeast')