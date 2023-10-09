clc
clear
close all
%% colors for plotting ber and sinr plots
colors = ["r", "b", "g", "k"];
%% Defining SNR values in range [0 20dB]
SNR_min = 0;
SNR_step = 1;
SNR_max = 15;
SNR_db = SNR_min:SNR_step:SNR_max;
SNR = 10.^(SNR_db/10);
%% Defining x[n] and innovation filter
x = [1/12, 0.5, 5/6, 0.5, 1/12];
L = 2;
Xz_roots = roots(x.');

% hence X(z) = F(z)F(1/z)
% where F(z) = (1+z^-1)(1+0.26792*z^-1)
f = [1, 1-Xz_roots(end), -Xz_roots(end)];
f = sqrt(x(1) /f(end))*f;
%% generate bpsk data
modulation = "BPSK";

% Choose equalizer type: ZF , MMSE or DFE
equalizer_type = "MMSE";

M = 2;
Es = 1; % E{|I_n|^2}

% equalizer parameters
taps_set = [5, 9, 13];
taps_set_len = length(taps_set);

% defining ber_vals
ber_vals = zeros(taps_set_len, length(SNR));
SINR_vals = zeros(taps_set_len, length(SNR));
ber_vals_theory = zeros(taps_set_len, length(SNR));
SINR_vals_theory = zeros(taps_set_len, length(SNR));

% for loop for calculating BER and SINR for different SNRs
for num = 1:taps_set_len
ntap = taps_set(num);
if ntap == inf
    ntap = 51;
end
k = (ntap-1)/2;
if equalizer_type == "DFE"
    k1 = (ntap-1)/2;
    k2 = k1;
end
if equalizer_type == "ZF"
    d = Equalizer(equalizer_type, ntap, x, 0, Es, 0);
    fprintf("ZF equalizer coefficients for %d taps:\n", ntap);
    disp(d);
    d0 = d(k+1);
    q = conv(x, d);
    q0 = q((length(q)+1)/2);
    qlen = length(q);
    permuts = 2*(dec2bin(2^(qlen)-1:-1:0)-'0')-1;
end

% generate
N = 1e6;  % Number of all bpsk data
I = generate_data(modulation, N);

for i = 1:length(SNR)
    % N0 = 1/(2*snr)
    snr = SNR(i);
    snr_sym = snr*log2(M);
    N0 = 1/(2*snr_sym);
    
    switch equalizer_type
        case "MMSE"
            d = Equalizer(equalizer_type, ntap, x, N0, Es, 0);
            d0 = d(k+1);
            q = conv(x, d);
            qlen = length(q);
            q0 = q((length(q)+1)/2);
            permuts = 2*(dec2bin(2^(qlen)-1:-1:0)-'0')-1;
        case "DFE"
            d = Equalizer(equalizer_type, ntap, x, N0, Es, [k1, k2]);
            df = d(1:(k1+1));
            d0 = df(end);
            db = d((k1+2):end);
            q = conv(x, df);
            qlen = length(q);
            q0 = q(end-L);
    end
    
    % generate noise
    noise = sqrt(N0)*(randn(1,N+L)+1i*randn(1,N+L));
    noise_out = conv(noise, f);
    nu = noise_out((L+1):(end-L));
    y = conv(x, I);
    y = y((L+1):(end-L)) + nu;
    
    % equalizer output
    switch equalizer_type
        case {"MMSE", "ZF"}
            Ihat = conv(d, y);
            Ihat = Ihat((k+1):(end-k));
            estimated_I = 2*(real(Ihat) > 0)-1;
            
            % calculate SINR
            SINR_num = sum(q0^2*(abs(I).^2))/N;
            SINR_den = sum(abs(Ihat-q0*I).^2)/N;
            SINR_vals(num, i) = SINR_num / SINR_den;
            SINR_vals_theory(num, i) = q0/(1-q0);
            
            % calculate simulation ber
            ber_vals(num, i) = sum(estimated_I ~= I)/N;
            
            % calculate ber in theory
            vec = conv(q, flip(d));
            w = vec((length(vec)+1)/2);
            ber_vals_theory(num, i) = (1/(2^(qlen)))*sum(qfunc(sqrt(Es*(permuts*q.').^2/(N0*w))));

            %ber_vals_theory(num, i) = qfunc(sqrt(SINR_vals(num, i)));
            
        case "DFE"
            SINR_num = 0;
            SINR_den = 0;
            estimated_I = zeros(1, N);
            idx = 1;
            for n = (k1+1):N
                idx_min_f = idx;
                idx_max_f = min(N, idx+k1);
                index_f = idx_min_f : idx_max_f;
                idx_min_b = max(1, idx-k2);
                idx_max_b = max(1, idx-1);
                index_b = idx_min_b : idx_max_b;

                yvec_f = zeros(1, k1+1);
                Ivec_b = zeros(1, k2);

                yvec_f((idx_min_f-idx+1):(idx_max_f-idx+1)) = y(idx_min_f:idx_max_f);
                Ivec_b((idx_min_b-idx+k2+1):(idx_max_b-idx+k2+1)) = estimated_I(idx_min_b:idx_max_b);
                if idx == N+50
                    dumb = 50;
                end
                if idx == 1
                    equ_out = fliplr(yvec_f)*df.';
                else
                    equ_out = fliplr(yvec_f)*df.' + fliplr(Ivec_b)*db.';
                end
                SINR_num = SINR_num + (1/N)*q0^2*abs(I(idx))^2;
                SINR_den = SINR_den + (1/N)*abs(equ_out - q0*I(idx))^2;
                estimated_I(idx) = (real(equ_out) > 0)*2-1;
                ber_vals(num, i) = ber_vals(num, i) + (I(idx) ~= estimated_I(idx))/N;
                idx = idx + 1;
            end
            
            SINR_vals(num, i) = SINR_num / SINR_den;
            SINR_vals_theory(num, i) = q0/(1-q0);
            ber_vals_theory(num, i) = qfunc(sqrt(SINR_vals(num, i)));
    end
end
end
%% Plotting BER
figure
for num = 1:taps_set_len
    semilogy(SNR_db, ber_vals(num,:), strcat(colors(num),'o-'), "DisplayName", strcat("Simulation and exact BER with ", num2str(taps_set(num)), " taps"));
    hold on
    % plot ber in theory
    semilogy(SNR_db, ber_vals_theory(num,:), strcat(colors(num),'s--'), "DisplayName", strcat("Theory BER with ", num2str(taps_set(num)), " taps"));
    grid on
    legend show
end
title(strcat("BPSK BER using ", equalizer_type ," equalization in range [0 15dB]"));
xlabel("$\frac{E\{|I_n|^2\}}{2N_0}$", "interpreter", "latex");
ylabel("Pe");
hold off
%% Plotting SINR
figure
for num = 1:taps_set_len
    semilogy(SNR_db, SINR_vals(num,:), strcat(colors(num),'o-'), "DisplayName", strcat("Simulation and exact SINR with ", num2str(taps_set(num)), " taps"));
    hold on
    % plot SINR in theory
    semilogy(SNR_db, SINR_vals_theory(num,:), strcat(colors(num),'s--'), "DisplayName", strcat("Theory SINR with ", num2str(taps_set(num)), " taps"));
    grid on
    legend show
end
title(strcat("BPSK SINR using ", equalizer_type ," equalization in range [0 15dB]"));
xlabel("$\frac{E\{|I_n|^2\}}{2N_0}$", "interpreter", "latex");
ylabel("SINR");
hold off