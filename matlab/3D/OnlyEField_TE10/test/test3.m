clear, clf
c0 = 299792458;
f = (100:2000)*10^6;
w = f*2*pi;
k0 = w*(1/c0);

a = 0.2;
gamma1 = 1j*sqrt(k0.^2-(pi/a).^2);

gamma2 = sqrt(((pi/a).^2)-k0.^2);
figure(1)
plot(f*10^-6,real(gamma1),'DisplayName','gamma1 real')
hold on
plot(f*10^-6,imag(gamma1),'DisplayName','gamma1 imag')
plot(f*10^-6,real(gamma2),'DisplayName','gamma2 real')
plot(f*10^-6,imag(gamma2),'DisplayName','gamma2 imag')
legend()
figure(2)

plot(f,k0)
%plot([(pi/a).^2,(pi/a).^2],[0,max(imag(gamma1))])
