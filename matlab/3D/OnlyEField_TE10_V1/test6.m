e10 = 0.0127;

f = (0.6:0.05:1.7)*10^9;
w = f*2*pi;
k0 = w.*(1/c0);
k_z10 = sqrt((k0.^2)-(pi/a).^2);
gamma = 1j*k_z10;

figure(1), clf, hold on
plot(f,real(k_z10*e10*2))
plot(f,imag(k_z10*e10*2))