set(0, 'DefaultTextInterpreter', 'none')
set(0, 'DefaultLegendInterpreter', 'none')
set(0, 'DefaultAxesTickLabelInterpreter', 'none')

figure(1),clf
subplot(2,1,1), hold on
plot(f,10*log10(DIO),'DisplayName','S11')
plot(f,10*log10(STON),'DisplayName','S12')
xlabel('Freqency [GHz]')
ylabel('dB')

lg = legend;
lg.FontSize = 11.5;
lg.Location = 'northwest';


subplot(2,1,2), hold on
plot(f,ST,'DisplayName','S11')
plot(f,VarN,'DisplayName','S12')

xlabel('Freqency [GHz]')
ylabel('Phase [$\degree$]')
set(gcf,'Position',[100 100 600 300])
saveas(gcf,'Sparameters_CST_FEM','svg')