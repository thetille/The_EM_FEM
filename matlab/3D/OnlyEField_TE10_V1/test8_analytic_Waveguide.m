clear
figure(1)
figure_offset = 0;
S21 = [];
for f = (0.1:0.005:1.5)*10^9
    if f > 0.75*10^9, figure_offset = 1; end
    res = 0.01;
    m = 1;
    n = 0;
    omega = f*2*pi;
    mu = 4*pi*10^(-7);
    eps = 8.8541878128*10^(-12);
    b = 0.1;
    a = 0.2;
    kc = sqrt(((m*pi/a)^2)+((n*pi/b)^2));
    k = omega*sqrt(mu*eps);
    be = sqrt((k^2) - (kc^2));
    
%     eta = sqrt(mu/eps);
%     Z = ((k*eta)/be);
%     E0 = sqrt(0.5*Z);
    
    
    A = (sqrt(2)*pi)/((a^(2/3))*sqrt(b)*sqrt(mu)*sqrt(omega)*sqrt(be));
    
    %(omega*mu*(a^3)*(A^2)*b)/(
    
    x = 0.1;
    y = 0.05;
    z = 0:res:0.4;    
    Ex_ana = 1i*(( 1j*omega*mu*n*pi)/((kc^2)*b)) *A* cos((m*pi*x)/a) * sin((n*pi*y)/b) * exp(1j*be.*z);
    Ey_ana = 1i*((-1j*omega*mu*m*pi)/((kc^2)*a)) *A* sin((m*pi*x)/a) * cos((n*pi*y)/b) * exp(1j*be.*z);
    Ez_ana = zeros(size(Ey_ana));
    
%     figure(1+figure_offset)
%     subplot(2,1,1), hold on
%     plot(z,abs(Ey_ana))
%     subplot(2,1,2), hold on
%     plot(z,angle(Ey_ana))
    
    
    S21 = [S21, (A*mu*omega*1i*exp(0.4*be*1i)*(0.1*cos(pi)-0.1))/kc^2];
    
end
figure(3)
plot((0.1:0.005:1.5),10*log10(abs(S21)))