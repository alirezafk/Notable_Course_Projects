clc
close all
N=10^5;     %Number of inputs of channel
X=rand(1,N);    %producing N standard uniform random variables
s=zeros(1,N);   %s is the input vector
n=randn(1,N);   %producing N standard normal independent noises
y_dB=0:0.1:18;  %The range of the plots
y=10.^(y_dB/10);    
sym_error=zeros(1,length(y));   %the vector of symbol error
bit_error=zeros(1,length(y));   %the vector of bit error

for M=[2 4 8 16]    %This for plots the probability error plots for different values of M. note that this code can be run for any M=2^k
    
Gray=grayCodes(log2(M));    %Gray is M*(log2(M)) matrix which ith row corresponds to the ith gray code, note that grayCodes is a function which gives a input n and the output is the gray codes of length n
for i=1:N                   %This for produces N random variables with M equiprobable outputs
    for j=0:(M-1)
     if((j/M<=X(i))&&(X(i)<=(j+1)/M))
            s(i)=2*j-(M-1);
     end
    end
end

C=sqrt((6*log2(M)/(M^2-1)).*y);     %The coefficients of of Am (note that we have normalized the equation r=s+n by multiplying 1/sqrt(2*N0), so n becomes a standard normal

for j=1:length(y)   %This for obtains the sym_error and bit_error vectors
    for i=1:N
        r=C(j)*s(i)+n(i);
        if(r>=(M-2)*C(j))
            if(s(i)~=(M-1))
                sym_error(j)=sym_error(j)+1/N;
                bit_error(j)=bit_error(j)+sum(xor(Gray(M,:),Gray((s(i)+M+1)/2,:)))/(N*log2(M));
            end
        elseif(r<=-(M-2)*C(j))
            if(s(i)~=-(M-1))
                sym_error(j)=sym_error(j)+1/N;
                bit_error(j)=bit_error(j)+sum(xor(Gray(1,:),Gray((s(i)+M+1)/2,:)))/(N*log2(M));
            end
        else
            z=floor(r/C(j));
            if(mod(z,2)==1)
                if(s(i)~=z)
                    sym_error(j)=sym_error(j)+1/N;
                    bit_error(j)=bit_error(j)+sum(xor(Gray((z+M+1)/2,:),Gray((s(i)+M+1)/2,:)))/(N*log2(M));
                end
            else
                if(s(i)~=(z+1))
                    sym_error(j)=sym_error(j)+1/N;
                    bit_error(j)=bit_error(j)+sum(xor(Gray((z+M+2)/2,:),Gray((s(i)+M+1)/2,:)))/(N*log2(M));
                end
            end
        end
    end
end

Analytical_symbol_Error=2*((M-1)/M)*qfunc(sqrt((6.*log2(M)/(M^2-1)).*y));
Analytical_bit_Error=Analytical_symbol_Error/log2(M);

figure
semilogy(y_dB,sym_error,'b','displayname',['Simulation symbol error M=' num2str(M)]);
axis([0 16 10^(-6) 1]);
hold on
semilogy(y_dB,Analytical_symbol_Error,'r','displayname',['analytical symbol error M=' num2str(M)]);
axis([0 18 10^(-6) 1]);
legend show
title(['M-PAM Symbol error for M=' num2str(M)]);

figure
semilogy(y_dB,bit_error,'b','displayname',['Simulation bit error M=' num2str(M)]);
axis([0 18 10^(-6) 1]);
hold on
semilogy(y_dB,Analytical_bit_Error,'r','displayname',['analytical bit error M=' num2str(M)]);
axis([0 18 10^(-6) 1]);
legend show
title(['M-PAM bit error for M=' num2str(M)]);
end
