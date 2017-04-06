
%                0.1nm lig  0.5nm lig   1nm lig   4nm lig 
dustdep_frac = [   0.01,       0.01,     0.01,     0.01;...
                   0.10,       0.10,     0.10,     0.10;...
                   0.25,       0.25,     0.25,     0.25;...
                   0.50,       0.50,     0.50,     0.50;...
                   1.00,       1.00,     1.00,     1.00;...
                   2.00,       2.00,     2.00,     2.00;...
                   4.00,       4.00,     4.00,     4.00;...
                   5.00,       5.00,     5.00,     5.00;...
                   10.0,       10.0,     10.0,     10.0];
% Gmol/yr               
dustdep_mass = [   0.033209,   0.033209, 0.033209, 0.033209;...
                   0.33209 ,   0.33209 , 0.33209 , 0.33209 ;...
                   0.83021 ,   0.83021 , 0.83021 , 0.83021 ;...
                   1.6604  ,   1.6604  , 1.6604  , 1.6604  ;...
                   3.3209  ,   3.3209  , 3.3209  , 3.3209  ;...
                   6.6417  ,   6.6417  , 6.6417  , 6.6417  ;...
                   13.283  ,   13.283  , 13.283  , 13.283  ;...
                   16.604  ,   16.604  , 16.604  , 16.604  ;...
                   33.209  ,   33.209  , 33.209  , 33.209  ];
                
supl_frac   = [    3.9279  ,   3.9279  , 3.9279  , 3.9279  ;... % somblgmt
                   0.389   ,   0.389   , 0.389   , 0.389   ];   % sombdc
% Gmol/yr               
supl_mass   = [    13.044  ,   13.044  , 13.044  , 13.044  ;... % somblgmt
                   1.2917  ,   1.2917  , 1.2917  , 1.2917  ];   % sombdc

%                0.1nm lig      0.5nm lig       1nm lig     4nm lig 
dustdep_pco2i = [2.7860513774131675e-04,2.8395671452335062e-04,2.7827560336940959e-04,2.7804760178173178e-04;...
                 2.7860513774131675e-04,2.8395671452335062e-04,2.7827560336940959e-04,2.7804760178173178e-04;...
                 2.7860513774131675e-04,2.8395671452335062e-04,2.7827560336940959e-04,2.7804760178173178e-04;...
                 2.7860513774131675e-04,2.8395671452335062e-04,2.7827560336940959e-04,2.7804760178173178e-04;...
                 2.7827560336940959e-04,2.7827560336940959e-04,2.7827560336940959e-04,2.7827560336940959e-04;...
                 2.7860513774131675e-04,2.8395671452335062e-04,2.7827560336940959e-04,2.7804760178173178e-04;...
                 2.7860513774131675e-04,2.8395671452335062e-04,2.7827560336940959e-04,2.7804760178173178e-04;...
                 2.7860513774131675e-04,2.8395671452335062e-04,2.7827560336940959e-04,2.7804760178173178e-04;...
                 2.7860513774131675e-04,2.8395671452335062e-04,2.7827560336940959e-04,2.7804760178173178e-04];

dustdep_pco2e = [2.7865692521628744e-04,2.8737397488010194e-04,3.1602158553843934e-04,3.0391382206044007e-04;...
                 2.7865070785807459e-04,2.8680560330299849e-04,3.0237666161925968e-04,2.7803668025404266e-04;...
                 2.7864052600680688e-04,2.8596643595837142e-04,2.8379960728898217e-04,2.7803668025404266e-04;...
                 2.7862407344460217e-04,2.8505872320980783e-04,2.7982528361834261e-04,2.7803668025404266e-04;...
                 2.7860513774131675e-04,2.8395671452335062e-04,2.7827560336940959e-04,2.7804760178173178e-04;...
                 2.7851454017573232e-04,2.8243927127074870e-04,2.7806430531225439e-04,2.7803668025404266e-04;...
                 2.7835124001088593e-04,2.8072225977773037e-04,2.7805307078012390e-04,2.7803668025404266e-04;...
                 2.7834765266869770e-04,2.8053255984440429e-04,2.7805083102920585e-04,2.7803668025404266e-04;...
                 2.7834194770564430e-04,2.8015863208790971e-04,2.7804759921983664e-04,2.7803668025404266e-04];

%                    
supl_pco2i = [      2.7860513774131675e-04,	2.8395671452335062e-04,	2.7827560336940959e-04,	2.7804760178173178e-04;... % somblgmt
                    2.7860513774131675e-04, 2.8395671452335062e-04,	2.7827560336940959e-04,	2.7804760178173178e-04];     % sombdc

supl_pco2e = [      2.7834416771106371e-04, 2.7804990487806218e-04,	2.7804812048034909e-04,	2.7803668025404266e-04;... % somblgmt
                    2.7864320877865286e-04, 2.8478539121352595e-04,	2.8027693996158463e-04,	2.7803668025404266e-04];     % sombdc


pco2_anomaly=dustdep_pco2e-dustdep_pco2i;
pco2_anomaly(5,:)=0; % Zero outthe control
supl_pco2_anomaly=supl_pco2e-supl_pco2i;

C = linspecer(4); % Colormap for the lines

figure
plot(dustdep_frac(:,2),dustdep_pco2e(:,2).*1e6,'o-','color',C(1,:),'MarkerEdgeColor',C(1,:),'MarkerFaceColor',C(1,:),'LineWidth',4)
hold on
plot(dustdep_frac(:,3),dustdep_pco2e(:,3).*1e6,'o-','color',C(3,:),'MarkerEdgeColor',C(3,:),'MarkerFaceColor',C(3,:),'LineWidth',4)
plot(dustdep_frac(:,4),dustdep_pco2e(:,4).*1e6,'o-','color',C(4,:),'MarkerEdgeColor',C(4,:),'MarkerFaceColor',C(4,:),'LineWidth',4)
% Add the spatially varying experiments
scatter(supl_frac(:,2),supl_pco2e(:,2).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(1,:),'LineWidth',1.5)
scatter(supl_frac(:,3),supl_pco2e(:,3).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(3,:),'LineWidth',1.5)
scatter(supl_frac(:,4),supl_pco2e(:,4).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(4,:),'LineWidth',1.5)
legend('0.5nm Ligand','1nm Ligand','4nm Ligand')
set(gca,'FontSize',16)
xlabel('Ratio of total dust deposition relative to preindustrial','FontSize',16)
ylabel('Atmospheric pCO_2 [uatm]','FontSize',16)
orient landscape
print -dpsc2 ~/Dropbox_Work/PIRES/pires_dust_only.ps

figure
plot(dustdep_frac(:,2),pco2_anomaly(:,2).*1e6,'o-','color',C(1,:),'MarkerEdgeColor',C(1,:),'MarkerFaceColor',C(1,:),'LineWidth',4)
hold on
plot(dustdep_frac(:,3),pco2_anomaly(:,3).*1e6,'o-','color',C(3,:),'MarkerEdgeColor',C(3,:),'MarkerFaceColor',C(3,:),'LineWidth',4)
plot(dustdep_frac(:,4),pco2_anomaly(:,4).*1e6,'o-','color',C(4,:),'MarkerEdgeColor',C(4,:),'MarkerFaceColor',C(4,:),'LineWidth',4)
% Add the spatially varying experiments
scatter(supl_frac(:,2),supl_pco2_anomaly(:,2).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(1,:),'LineWidth',1.5)
scatter(supl_frac(:,3),supl_pco2_anomaly(:,3).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(3,:),'LineWidth',1.5)
scatter(supl_frac(:,4),supl_pco2_anomaly(:,4).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(4,:),'LineWidth',1.5)
legend('0.5nm Ligand','1nm Ligand','4nm Ligand')
set(gca,'FontSize',16)
xlabel('Ratio of total dust deposition relative to preindustrial','FontSize',16)
ylabel('Atmospheric pCO_2 anomaly [uatm]','FontSize',16)
orient landscape
print -dpsc2 ~/Dropbox_Work/PIRES/pires_dust_only_anomaly.ps

figure
plot(dustdep_mass(:,2),pco2_anomaly(:,2).*1e6,'o-','color',C(1,:),'MarkerEdgeColor',C(1,:),'MarkerFaceColor',C(1,:),'LineWidth',4)
hold on
plot(dustdep_mass(:,3),pco2_anomaly(:,3).*1e6,'o-','color',C(3,:),'MarkerEdgeColor',C(3,:),'MarkerFaceColor',C(3,:),'LineWidth',4)
plot(dustdep_mass(:,4),pco2_anomaly(:,4).*1e6,'o-','color',C(4,:),'MarkerEdgeColor',C(4,:),'MarkerFaceColor',C(4,:),'LineWidth',4)
% Add the spatially varying experiments
scatter(supl_mass(:,2),supl_pco2_anomaly(:,2).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(1,:),'LineWidth',1.5)
scatter(supl_mass(:,3),supl_pco2_anomaly(:,3).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(3,:),'LineWidth',1.5)
scatter(supl_mass(:,4),supl_pco2_anomaly(:,4).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(4,:),'LineWidth',1.5)
legend('0.5nm Ligand','1nm Ligand','4nm Ligand')
set(gca,'FontSize',16)
xlabel('Iron/Dust deposition [Gmol/yr]','FontSize',16)
ylabel('Atmospheric pCO_2 anomaly [uatm]','FontSize',16)
orient landscape
print -dpsc2 ~/Dropbox_Work/PIRES/pires_dust_only_mass_anomaly.ps


   


















