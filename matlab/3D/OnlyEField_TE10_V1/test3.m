clear

a = 0.2;
b = 0.1;
c0 = 299792458; % speed of light in vacuum

f = (0.6:0.005:1.7)*10^9;
w = f*2*pi;
k0 = w*(1/c0);
k_z10 = sqrt((k0.^2)-(pi/a).^2);

figure(1), clf
plot(f,real(k_z10))
hold on
plot(f,imag(k_z10))
plot(f,k0)

length = 0.4;
figure(2), clf
plot(f,real(k_z10).*length./(2*pi))
hold on
plot(f,imag(k_z10).*length./(2*pi))
plot(f,k0.*length./(2*pi))

k_z10f = @(f) sqrt((((f.*2.*pi)*(1./299792458)).^2)-(pi./0.2).^2).*length./(2*pi);

%k_z10f([0.78416,0.88,1.02,1.185,1.370,1.565]*10^9)
k_z10f([0.75,0.837984,1.055,1.345,1.665,]*10^9)