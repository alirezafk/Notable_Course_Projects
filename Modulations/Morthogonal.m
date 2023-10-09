clc
close all
N=10^4;     %Number of inputs of channel
X=rand(1,N);        %producing N standard uniform random variables
y_dB=0:0.1:18;      %The range of the plots
y=10.^(y_dB/10);
Analytical_symbol_Error=zeros(1,length(y));
Analytical_bit_Error=zeros(1,length(y));
graycode=zeros(1,N);
sym_error=zeros(1,length(y));       %the vector of symbol error
bit_error=zeros(1,length(y));       %the vector of bit error

for M=[4 8 16]       %This for plots the probability error plots for different values of M. note that this code can be run for any M=2^k

Gray=grayCodes(log2(M));    %Gray is M*(log2(M)) matrix which ith row corresponds to the ith gray code, note that grayCodes is a function which gives a input n and the output is the gray codes of length n
s=zeros(M,N);       %s is the input vector
r=zeros(M,1);
n=randn(M,N);       %producing N standard normal independent noises
d=zeros(1,M);
signals=eye(M);

for i=1:N       %This for produces N random variables with M equiprobable outputs
    for j=0:(M-1)
     if((j/M<=X(i))&&(X(i)<=(j+1)/M))
        s(:,i)=signals(:,j+1);
        graycode(i)=j+1;
     end
    end
end

C=sqrt(2*log2(M).*y);       %The coefficients of of Am (note that we have normalized the equation r=s+n by multiplying 1/sqrt(2*N0), so n becomes a standard normal

for j=1:length(y)
    for i=1:N
        r=C(j)*s(:,i)+n(:,i);
        for t=1:M
            d(t)=norm(r-C(j)*signals(:,t));
        end
        Index=find(d==min(d));
        correct=0;
        for v=1:length(Index)
            if(signals(:,Index(v))==s(:,i))
                correct=1;
            end
        end
        if(correct==0)
            sym_error(j)=sym_error(j)+1/N;
            bit_error(j)=bit_error(j)+sum(xor(Gray(graycode(i)),Gray(Index(1))))/(N*log2(M));
        end
    end
end

for i=1:length(y)
    I=@(h)(1-(1-qfunc(sqrt(2*log2(M)*y(i))+h)).^(M-1)).*(1/sqrt(2*pi)).*exp(-(h.^2)/2);
    Analytical_symbol_Error(i)=integral(I,-20,20);
end

Analytical_bit_Error=Analytical_symbol_Error/log2(M);
figure
semilogy(y_dB,sym_error,'b','displayname',['Simulation symbol error M=' num2str(M)]);
axis([0 16 10^(-8) 1]);
hold on
semilogy(y_dB,Analytical_symbol_Error,'r','displayname',['analytical symbol error M=' num2str(M)]);
axis([0 18 10^(-8) 1]);
legend show
title(['M-orthogonal Symbol error for M=' num2str(M)]);

figure
semilogy(y_dB,bit_error,'b','displayname',['Simulation bit error M=' num2str(M)]);
axis([0 18 10^(-8) 1]);
hold on
semilogy(y_dB,Analytical_bit_Error,'r','displayname',['analytical bit error M=' num2str(M)]);
axis([0 18 10^(-8) 1]);
legend show
title(['M-orthogonal bit error for M=' num2str(M)]);
end
