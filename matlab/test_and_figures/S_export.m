clear

set(0, 'DefaultTextInterpreter', 'none')
set(0, 'DefaultLegendInterpreter', 'none')
set(0, 'DefaultAxesTickLabelInterpreter', 'none')

%%%%% warning wrong labeles/names in data CST %%%%%%%
cols = get(groot,'DefaultAxesColorOrder');

CST_FEM_color = cols(1,:);
CST_FIT_color = cols(2,:);
Our_FEM_color = cols(3,:);

CST_FEM = load('waveguide_with_ports');
CST_FIT = load('waveguide_with_ports_FIT');
Our_FEM = load('../3D/OnlyEField_TE10_V1/res/waveguide_model3_flat6/test2/Sparamters');

S21 = [];
ana_f = (0.1:0.005:1.5)*10^9;
for f = ana_f
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
    A = (sqrt(2)*pi)/((a^(2/3))*sqrt(b)*sqrt(mu)*sqrt(omega)*sqrt(be));
    S21 = [S21, (A*mu*omega*1i*exp(0.4*be*1i)*(0.1*cos(pi)-0.1))/kc^2];
end



figure(1),clf
t = tiledlayout(2,1);
ax = nexttile; 
hold on
plot(CST_FEM.f,10*log10(CST_FEM.S12_Mag),'DisplayName','CST FEM S11','color',CST_FEM_color)
plot(CST_FIT.f,10*log10(CST_FIT.S11_Mag),'DisplayName','CST FIT S11','color',CST_FIT_color)
plot(Our_FEM.f_list*10^-9,10*log10(abs(Our_FEM.S_par(:,1,1))),'DisplayName','Our FEM S11','color',Our_FEM_color)
plot([0.74948, 0.74948],[-80,0],'--','color',[0.5,0.5,0.5],'DisplayName','Cut off frequency')

set(ax, 'XTick', 0:0.15:1.5)

grid on
xlabel('Freqency [GHz]')
ylabel('dB')


lh =legend(ax,'Location','NorthOutside','Orientation','Horizontal');
lh.FontSize = 11.5;



ax = nexttile; hold on
plot(CST_FEM.f,CST_FEM.S11_Phase,'DisplayName','CST FEM S11','color',CST_FEM_color)
plot(CST_FIT.f,CST_FIT.S11_Phase,'DisplayName','CST FIT S11','color',CST_FIT_color)
plot(Our_FEM.f_list*10^-9,angle(Our_FEM.S_par(:,1,1))*180/pi,'DisplayName','Our FEM S11','color',Our_FEM_color)
xlabel('Freqency [GHz]')
ylabel('Phase [$\degree$]')
plot([0.74948, 0.74948],[-180,180],'--','color',[0.5,0.5,0.5],'DisplayName','Cut off frequency')
set(ax, 'XTick', 0:0.15:1.5)
set(gcf,'Position',[100 100 700 460])
saveas(gcf,'Sparameters_compear_S11','svg')
%nexttile


figure(2),clf
t = tiledlayout(2,1);
ax = nexttile; hold on
plot(CST_FEM.f,10*log10(CST_FEM.S11_Mag),'DisplayName','CST FEM S21','color',CST_FEM_color)
plot(CST_FIT.f,10*log10(CST_FIT.S12_Mag),'DisplayName','CST FIT S21','color',CST_FIT_color)
plot(Our_FEM.f_list*10^-9,10*log10(abs(Our_FEM.S_par(:,2,1))),'DisplayName','Our FEM S21','color',Our_FEM_color)
plot(ana_f*10^-9,10*log10(abs(S21)),'DisplayName','Analytical S21')
plot([0.74948, 0.74948],[-80,0],'--','color',[0.5,0.5,0.5],'DisplayName','Cut off frequency')
set(ax, 'XTick', 0:0.15:1.5)
xlabel('Freqency [GHz]')
ylabel('dB')
grid on

lh =legend(ax,'Location','NorthOutside','Orientation','Horizontal');
lh.FontSize = 11.5;

ax = nexttile; hold on
plot(CST_FIT.f,CST_FIT.S12_phase,'DisplayName','CST FIT S12','color',CST_FIT_color)
plot(CST_FEM.f,CST_FEM.S12_Phase,'DisplayName','CST FEM S12','color',CST_FEM_color)
plot(Our_FEM.f_list*10^-9,angle(Our_FEM.S_par(:,2,1))*180/pi,'DisplayName','Our FEM S21','color',Our_FEM_color)
%plot(ana_f*10^-9,angle(S21),'DisplayName','Analytical S21')
set(ax, 'XTick', 0:0.15:1.5)
plot([0.74948, 0.74948],[-180,180],'--','color',[0.5,0.5,0.5],'DisplayName','Cut off frequency')
xlabel('Freqency [GHz]')
ylabel('Phase [$\degree$]')
set(gcf,'Position',[100 100 700 460])
saveas(gcf,'Sparameters_compear_S12','svg')


figure(3),clf,hold on
start = 464;
start2 = 1;
plot(CST_FEM.f(start:end),unwrap(CST_FEM.S12_Phase(start:end)*pi/180)*180/pi,'DisplayName','CST FEM S21','color',CST_FEM_color)
plot(CST_FIT.f(start:end),unwrap(CST_FIT.S12_phase(start:end)*pi/180)*180/pi,'DisplayName','CST FIT S21','color',CST_FIT_color)
plot(Our_FEM.f_list(start2:end)*10^-9,unwrap(angle(Our_FEM.S_par(start2:end,2,1)))*180/pi,'DisplayName','Our FEM S21 True','color',[0.8,0.8,0.8])
plot(Our_FEM.f_list(start2:end)*10^-9,-1*unwrap(angle(Our_FEM.S_par(start2:end,2,1)))*180/pi,'DisplayName','Our FEM S21 Flipped','color',Our_FEM_color)
plot([0.74948, 0.74948],[-180,180],'--','color',[0.5,0.5,0.5],'DisplayName','Cut off frequency')
lh = legend();
lh.FontSize = 11.5;
xlabel('Freqency [GHz]')
ylabel('Phase [$\degree$]')
xlim([0.74948,1.5]);
grid on
set(gcf,'Position',[100 100 700 460])
saveas(gcf,'Sparameters_compear_S21_unwrap','svg')

start = 460;
start3 = 2330;
start2 = 27;
figure(4), clf, hold on
plot(CST_FEM.f(start3:end),10*log10(CST_FEM.S11_Mag(start3:end)),'DisplayName','CST FEM S21','color',CST_FEM_color)
plot(CST_FIT.f(start:end),10*log10(CST_FIT.S12_Mag(start:end)),'DisplayName','CST FIT S21','color',CST_FIT_color)
plot(Our_FEM.f_list(start2:end)*10^-9,10*log10(abs(Our_FEM.S_par((start2:end),2,1))),'DisplayName','Our FEM S21','color',Our_FEM_color)
%plot(ana_f*10^-9,10*log10(abs(S21)),'DisplayName','Analytical S21')
%plot([0.74948, 0.74948],[-80,0],'--','color',[0.5,0.5,0.5],'DisplayName','Cut off frequency')
set(ax, 'XTick', 0:0.15:1.5)
xlabel('Freqency [GHz]')
ylabel('dB')
grid on
ylim([-0.015,0.015])
xlim([0.74948,1.5])
lh =legend(ax,'Location','NorthOutside','Orientation','Horizontal');
lh.FontSize = 11.5;
legend()

out_FEM_new = interp1(Our_FEM.f_list*10^-9,abs(Our_FEM.S_par(:,2,1)),CST_FIT.f);
CST_FEM_new = interp1(CST_FEM.f,CST_FEM.S11_Mag,CST_FIT.f);
pass_at = 465;

mean(out_FEM_new-CST_FEM_new);
%fprinf('Erorr S12 in pass:' max(CST_FEM.S11_Mag-10*log10(abs(Our_FEM.S_par(:,2,1)))))