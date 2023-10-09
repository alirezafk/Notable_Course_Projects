clc
close all
N=10^5;     %Number of inputs of channel
X=rand(1,N);        %producing N standard uniform random variables
s=zeros(2,N);       %s is the input vector
r=zeros(2,1);
u=zeros(2,1);
n=randn(2,N);       %producing N standard normal independent noises
graycode=zeros(1,N);
y_dB=0:0.1:18;      %The range of the plots
y=10.^(y_dB/10);
sym_error=zeros(1,length(y));       %the vector of symbol error
bit_error=zeros(1,length(y));       %the vector of bit error

for M=[4 16]         %This for plots the probability error plots for different values of M. note that this code can be run for any M=2^k

Gray=grayCodes(log2(M));    %Gray is M*(log2(M)) matrix which ith row corresponds to the ith gray code, note that grayCodes is a function which gives a input n and the output is the gray codes of length n
RM=sqrt(M);

for i=1:N       %This for produces N random variables with M equiprobable outputs
    for j=0:(M-1)
     if((j/M<=X(i))&&(X(i)<=(j+1)/M))
         s(1,i)=-(RM-1)+2*mod(j,RM);
         s(2,i)=-(RM-1)+2*floor(j/RM);
         graycode(i)=j+1;
     end
    end
end

C=sqrt((3*log2(M)/(M-1)).*y);       %The coefficients of of Am (note that we have normalized the equation r=s+n by multiplying 1/sqrt(2*N0), so n becomes a standard normal

for j=1:length(y)
    for i=1:N
        r=C(j)*s(:,i)+n(:,i);
        z1=floor(r(1,1)/C(j));
        z2=floor(r(2,1)/C(j));
        if(mod(z1,2)==0)
            z1=z1+1;
        end
        if(mod(z2,2)==0)
            z2=z2+1;
        end
        if(z1>=RM-1)
            z1=RM-1;
        end
        if(z1<=-(RM-1))
            z1=-(RM-1);
        end
        if(z2>=RM-1)
            z2=RM-1;
        end
        if(z2<=-(RM-1))
            z2=-(RM-1);
        end
        gray_index=1+((z2+RM-1)/2)*RM+(z1+RM-1)/2;
        if(z1~=s(1,i)||z2~=s(2,i))
            sym_error(j)=sym_error(j)+1/N;
            bit_error(j)=bit_error(j)+sum(xor(Gray(graycode(i)),Gray(gray_index)))/(log2(M)*N);
        end    
    end
end

PE=2*((RM-1)/RM)*qfunc(sqrt((3.*log2(M)/(M-1)).*y));
Analytical_symbol_Error=2*PE-PE.^2;
Analytical_bit_Error=Analytical_symbol_Error/log2(M);

figure
semilogy(y_dB,sym_error,'b','displayname',['Simulation symbol error M=' num2str(M)]);
axis([0 16 10^(-8) 1]);
hold on
semilogy(y_dB,Analytical_symbol_Error,'r','displayname',['analytical symbol error M=' num2str(M)]);
axis([0 18 10^(-8) 1]);
legend show
title(['M-QAM Symbol error for M=' num2str(M)]);

figure
semilogy(y_dB,bit_error,'b','displayname',['Simulation bit error M=' num2str(M)]);
axis([0 18 10^(-8) 1]);
hold on
semilogy(y_dB,Analytical_bit_Error,'r','displayname',['analytical bit error M=' num2str(M)]);
axis([0 18 10^(-8) 1]);
legend show
title(['M-QAM bit error for M=' num2str(M)]);
end
