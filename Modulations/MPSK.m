clc
close all
N=10^5;     %Number of inputs of channel
X=rand(1,N);        %producing N standard uniform random variables
s=zeros(2,N);       %s is the input vector
s_theta=zeros(1,N);
r=zeros(2,1);
n=randn(2,N);       %producing N standard normal independent noises
y_dB=0:0.1:18;      %The range of the plots
y=10.^(y_dB/10);
sym_error=zeros(1,length(y));       %the vector of symbol error
bit_error=zeros(1,length(y));       %the vector of bit error

for M=[2 4 8 16]         %This for plots the probability error plots for different values of M. note that this code can be run for any M=2^k
 
Gray=grayCodes(log2(M));    %Gray is M*(log2(M)) matrix which ith row corresponds to the ith gray code, note that grayCodes is a function which gives a input n and the output is the gray codes of length n

for i=1:N       %This for produces N random variables with M equiprobable outputs
    for j=0:(M-1)
     if((j/M<=X(i))&&(X(i)<=(j+1)/M))
        s(1,i)=cos(2*pi*(j/M));
        s(2,i)=sin(2*pi*(j/M));
        s_theta(i)=j;
     end
    end
end

C=sqrt((2*log2(M)).*y);     %The coefficients of of Am (note that we have normalized the equation r=s+n by multiplying 1/sqrt(2*N0), so n becomes a standard normal

for j=1:length(y)
    for i=1:N
        r(1,1)=C(j)*s(1,i)+n(1,i);
        r(2,1)=C(j)*s(2,i)+n(2,i);
        [theta,rho]=cart2pol(r(1,1),r(2,1));
        if(theta<0)
            theta=theta+2*pi;
        end
        if((theta<=pi/M&&theta>=0)||(theta>=(2*pi-pi/M)&&theta<2*pi))
            if(s_theta(i)~=0)
                sym_error(j)=sym_error(j)+1/N;
                bit_error(j)=bit_error(j)+sum(xor(Gray(1,:),Gray(s_theta(i)+1,:)))/N;
            end
        else
            z=floor(theta*(M/pi));
            if(mod(z,2)==0)
                if(s_theta(i)~=z/2)
                    sym_error(j)=sym_error(j)+1/N;
                    bit_error(j)=bit_error(j)+sum(xor(Gray(z/2+1,:),Gray(s_theta(i)+1,:)))/(N*log2(M));
                end
            else
                if(s_theta(i)~=(z+1)/2)
                    sym_error(j)=sym_error(j)+1/N;
                    bit_error(j)=bit_error(j)+sum(xor(Gray((z+1)/2+1,:),Gray(s_theta(i)+1,:)))/(N*log2(M));
                end 
            end
        end
    end
end

Analytical_symbol_Error=2*qfunc(sqrt(2.*log2(M).*(sin(pi/M))^2.*y));
Analytical_bit_Error=Analytical_symbol_Error/log2(M);

figure
semilogy(y_dB,sym_error,'b','displayname',['Simulation symbol error M=' num2str(M)]);
axis([0 16 10^(-8) 1]);
hold on
semilogy(y_dB,Analytical_symbol_Error,'r','displayname',['analytical symbol error M=' num2str(M)]);
axis([0 18 10^(-8) 1]);
legend show
title(['M-PSK Symbol error for M=' num2str(M)]);

figure
semilogy(y_dB,bit_error,'b','displayname',['Simulation bit error M=' num2str(M)]);
axis([0 18 10^(-8) 1]);
hold on
semilogy(y_dB,Analytical_bit_Error,'r','displayname',['analytical bit error M=' num2str(M)]);
axis([0 18 10^(-8) 1]);
legend show
title(['M-PSK bit error for M=' num2str(M)]);
end