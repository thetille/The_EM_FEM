e10 = 0.0127;
clear
a = 0.2;
b = 0.1;
c0 = 299792458;

f = (0.6:0.005:1.7)*10^9;
w = f*2*pi;
k0 = w.*(1/c0);
k_z10 = sqrt((k0.^2)-(pi/a).^2);
gamma = 1j*k_z10;

Z = (w.*(pi*4*10^(-7))./k_z10);
E0 = sqrt(0.5.*Z);

figure(1), hold on
plot(f*10^-9,real(Z))
plot(f*10^-9,imag(Z))

% figure(1), clf, hold on
% plot(f,real(k_z10*e10*2))
% plot(f,imag(k_z10*e10*2))