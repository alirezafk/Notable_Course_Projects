m=8;
N=10^6;
X=rand(1,N);
s=zeros(2,N);
s_theta=zeros(1,N);
r=zeros(2,1);
n=randn(2,N);
y_dB=0:0.1:16;
y=10.^(y_dB/10);
Error=zeros(1,length(y));

for i=1:N
    for j=0:(M-1)
     if((j/M<=X(i))&&(X(i)<=(j+1)/M))
        s(1,i)=cos(2*pi*(j/M));
        s(2,i)=sin(2*pi*(j/M));
        s_theta(i)=2*pi*(j/M);
     end
    end
end

C=sqrt((2*log2(M)).*y);

for j=1:length(y)
    for i=1:N
        r(1,1)=C(j)*s(1,i)+n(1,i);
        r(2,1)=C(j)*s(2,i)+n(2,i);
        [theta,rho]=cart2pol(r(1,1),r(2,1));
        if(theta<0)
            theta=theta+2*pi;
        end
        if(s_theta(i)==0)
             if((theta>pi/M)&&(theta<(2*pi-pi/M)))
                 Error(j)=Error(j)+1/N;
             end
        else
            if(~((theta>=s_theta(i)-pi/M)&&(theta<=s_theta(i)+pi/M)))
                Error(j)=Error(j)+1/N;
            end
        end
    end
end

Analytical_Error=2*qfunc(sqrt(2.*log2(M).*(sin(pi/M))^2.*y));

figure
semilogy(y_dB,Error);
hold on
semilogy(y_dB,Analytical_Error);
