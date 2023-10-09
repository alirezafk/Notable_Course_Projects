w=0:0.01:pi;
H1=(exp(-1i*7.*w)/15).*((sin(15.*w/2)./sin(w/2))-(0.5*exp(-1i*pi/15).*sin(15.*w/2))/(sin((w-2*pi/15)/2))-(0.5*exp(1i*pi/15).*sin(15.*w/2))/(sin((w+2*pi/15)/2)));
H2=(exp(-1i*7.*w)/15).*((sin(15.*w/2)./sin(w/2))+(0.5*exp(-1i*pi/15).*sin(15.*w/2))/(sin((w-2*pi/15)/2))+(0.5*exp(1i*pi/15).*sin(15.*w/2))/(sin((w+2*pi/15)/2)));

figure
plot(w,abs(H1));

figure
plot(w,abs(H2));