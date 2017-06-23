% 
% %                0.1nm lig  0.5nm lig   1nm lig   4nm lig 
% dustdep_frac = [   0.01,       0.01,     0.01,     0.01;...
%                    0.10,       0.10,     0.10,     0.10;...
%                    0.25,       0.25,     0.25,     0.25;...
%                    0.50,       0.50,     0.50,     0.50;...
%                    1.00,       1.00,     1.00,     1.00;...
%                    2.00,       2.00,     2.00,     2.00;...
%                    4.00,       4.00,     4.00,     4.00;...
%                    5.00,       5.00,     5.00,     5.00;...
%                    10.0,       10.0,     10.0,     10.0];
% % Gmol/yr               
% dustdep_mass = [   0.033209,   0.033209, 0.033209, 0.033209;...
%                    0.33209 ,   0.33209 , 0.33209 , 0.33209 ;...
%                    0.83021 ,   0.83021 , 0.83021 , 0.83021 ;...
%                    1.6604  ,   1.6604  , 1.6604  , 1.6604  ;...
%                    3.3209  ,   3.3209  , 3.3209  , 3.3209  ;...
%                    6.6417  ,   6.6417  , 6.6417  , 6.6417  ;...
%                    13.283  ,   13.283  , 13.283  , 13.283  ;...
%                    16.604  ,   16.604  , 16.604  , 16.604  ;...
%                    33.209  ,   33.209  , 33.209  , 33.209  ];
%                 
% supl_frac   = [    3.9279  ,   3.9279  , 3.9279  , 3.9279  ;... % somblgmt
%                    0.389   ,   0.389   , 0.389   , 0.389   ];   % sombdc
% % Gmol/yr               
% supl_mass   = [    13.044  ,   13.044  , 13.044  , 13.044  ;... % somblgmt
%                    1.2917  ,   1.2917  , 1.2917  , 1.2917  ];   % sombdc
% 
% %                0.1nm lig      0.5nm lig       1nm lig     4nm lig 
% dustdep_pco2i = [2.7860513774131675e-04,2.8395671452335062e-04,2.7827560336940959e-04,2.7804760178173178e-04;...
%                  2.7860513774131675e-04,2.8395671452335062e-04,2.7827560336940959e-04,2.7804760178173178e-04;...
%                  2.7860513774131675e-04,2.8395671452335062e-04,2.7827560336940959e-04,2.7804760178173178e-04;...
%                  2.7860513774131675e-04,2.8395671452335062e-04,2.7827560336940959e-04,2.7804760178173178e-04;...
%                  2.7827560336940959e-04,2.7827560336940959e-04,2.7827560336940959e-04,2.7827560336940959e-04;...
%                  2.7860513774131675e-04,2.8395671452335062e-04,2.7827560336940959e-04,2.7804760178173178e-04;...
%                  2.7860513774131675e-04,2.8395671452335062e-04,2.7827560336940959e-04,2.7804760178173178e-04;...
%                  2.7860513774131675e-04,2.8395671452335062e-04,2.7827560336940959e-04,2.7804760178173178e-04;...
%                  2.7860513774131675e-04,2.8395671452335062e-04,2.7827560336940959e-04,2.7804760178173178e-04];
% 
% dustdep_pco2e = [2.7865692521628744e-04,2.8737397488010194e-04,3.1602158553843934e-04,3.0391382206044007e-04;...
%                  2.7865070785807459e-04,2.8680560330299849e-04,3.0237666161925968e-04,2.7803668025404266e-04;...
%                  2.7864052600680688e-04,2.8596643595837142e-04,2.8379960728898217e-04,2.7803668025404266e-04;...
%                  2.7862407344460217e-04,2.8505872320980783e-04,2.7982528361834261e-04,2.7803668025404266e-04;...
%                  2.7860513774131675e-04,2.8395671452335062e-04,2.7827560336940959e-04,2.7804760178173178e-04;...
%                  2.7851454017573232e-04,2.8243927127074870e-04,2.7806430531225439e-04,2.7803668025404266e-04;...
%                  2.7835124001088593e-04,2.8072225977773037e-04,2.7805307078012390e-04,2.7803668025404266e-04;...
%                  2.7834765266869770e-04,2.8053255984440429e-04,2.7805083102920585e-04,2.7803668025404266e-04;...
%                  2.7834194770564430e-04,2.8015863208790971e-04,2.7804759921983664e-04,2.7803668025404266e-04];
% 
% %                    
% supl_pco2i = [      2.7860513774131675e-04,	2.8395671452335062e-04,	2.7827560336940959e-04,	2.7804760178173178e-04;... % somblgmt
%                     2.7860513774131675e-04, 2.8395671452335062e-04,	2.7827560336940959e-04,	2.7804760178173178e-04];     % sombdc
% 
% supl_pco2e = [      2.7834416771106371e-04, 2.7804990487806218e-04,	2.7804812048034909e-04,	2.7803668025404266e-04;... % somblgmt
%                     2.7864320877865286e-04, 2.8478539121352595e-04,	2.8027693996158463e-04,	2.7803668025404266e-04];     % sombdc

%                 1nm lig     1nm lig   1nm lig     1nm lig   1nm lig   1nm lig
%                 a2 dust     a4 dust   a6 dust     a2 all    a4 all    a6 all
dustdep_frac = [   0.00,       0.00,     0.00,       0.00,     0.00,     0.00;...
                   0.01,       0.01,     0.01,       0.01,     0.01,     0.01;...
                   0.10,       0.10,     0.10,       0.10,     0.10,     0.10;...
                   0.25,       0.25,     0.25,       0.25,     0.25,     0.25;...
                   0.50,       0.50,     0.50,       0.50,     0.50,     0.50;...
                   1.00,       1.00,     1.00,       1.00,     1.00,     1.00;...
                   2.00,       2.00,     2.00,       2.00,     2.00,     2.00;...
                   4.00,       4.00,     4.00,       4.00,     4.00,     4.00;...
                   5.00,       5.00,     5.00,       5.00,     5.00,     5.00;...
                   10.0,       10.0,     10.0,       10.0,     10.0,     10.0];
% Gmol/yr               
dustdep_mass = [   0.000000,   0.000000, 0.000000,   0.000000, 0.000000, 0.000000;...
                   0.033209,   0.033209, 0.033209,   0.033209, 0.033209, 0.033209;...
                   0.33209 ,   0.33209 , 0.33209 ,   0.33209 , 0.33209 , 0.33209 ;...
                   0.83021 ,   0.83021 , 0.83021 ,   0.83021 , 0.83021 , 0.83021 ;...
                   1.6604  ,   1.6604  , 1.6604  ,   1.6604  , 1.6604  , 1.6604  ;...
                   3.3209  ,   3.3209  , 3.3209  ,   3.3209  , 3.3209  , 3.3209  ;...
                   6.6417  ,   6.6417  , 6.6417  ,   6.6417  , 6.6417  , 6.6417  ;...
                   13.283  ,   13.283  , 13.283  ,   13.283  , 13.283  , 13.283  ;...
                   16.604  ,   16.604  , 16.604  ,   16.604  , 16.604  , 16.604  ;...
                   33.209  ,   33.209  , 33.209  ,   33.209  , 33.209  , 33.209  ];
                
supl_frac   = [    3.9279  ,   3.9279  , 3.9279  ,   3.9279  , 3.9279  , 3.9279  ;... % somblgmt
                   0.389   ,   0.389   , 0.389   ,   0.389   , 0.389   , 0.389   ];   % sombdc
% Gmol/yr               
supl_mass   = [    13.044  ,   13.044  , 13.044  ,   13.044  , 13.044  , 13.044  ;... % somblgmt
                   1.2917  ,   1.2917  , 1.2917  ,   1.2917  , 1.2917  , 1.2917  ];   % sombdc

%                        1nm lig                1nm lig                 1nm lig                 1nm lig                 1nm lig                 1nm lig  
%                        a2 dust                a4 dust                 a6 dust                 a2 all                  a4 all                   a6 all
dustdep_pco2i = [2.7827560336940959e-04, 2.7827560336940959e-04, 2.7827560336940959e-04, 2.7805059701209667e-04, 2.7827560336940959e-04, 2.7827560336940959e-04;...
                 2.7827560336940959e-04, 2.7827560336940959e-04, 2.7827560336940959e-04, 2.7805059701209667e-04, 2.7827560336940959e-04, 2.7827560336940959e-04;...
                 2.7827560336940959e-04, 2.7827560336940959e-04, 2.7827560336940959e-04, 2.7805059701209667e-04, 2.7827560336940959e-04, 2.7827560336940959e-04;...
                 2.7827560336940959e-04, 2.7827560336940959e-04, 2.7827560336940959e-04, 2.7805059701209667e-04, 2.7827560336940959e-04, 2.7827560336940959e-04;...
                 2.7827560336940959e-04, 2.7827560336940959e-04, 2.7827560336940959e-04, 2.7805059701209667e-04, 2.7827560336940959e-04, 2.7827560336940959e-04;...
                 2.7827560336940959e-04, 2.7827560336940959e-04, 2.7827560336940959e-04, 2.7805059701209667e-04, 2.7827560336940959e-04, 2.7827560336940959e-04;...
                 2.7827560336940959e-04, 2.7827560336940959e-04, 2.7827560336940959e-04, 2.7805059701209667e-04, 2.7827560336940959e-04, 2.7827560336940959e-04;...
                 2.7827560336940959e-04, 2.7827560336940959e-04, 2.7827560336940959e-04, 2.7805059701209667e-04, 2.7827560336940959e-04, 2.7827560336940959e-04;...
                 2.7827560336940959e-04, 2.7827560336940959e-04, 2.7827560336940959e-04, 2.7805059701209667e-04, 2.7827560336940959e-04, 2.7827560336940959e-04;...
                 2.7827560336940959e-04, 2.7827560336940959e-04, 2.7827560336940959e-04, 2.7805059701209667e-04, 2.7827560336940959e-04, 2.7827560336940959e-04];

dustdep_pco2e = [3.1735176639198438e-04, 3.4056542087654746e-04, 3.4635531093347229e-04, 2.9860463184297470e-04, 2.9480437170197233e-04, 2.9407540226501739e-04;...
                 3.1602158553843934e-04, 3.3634911785843758e-04, 3.4376981536857480e-04, 2.9703380236476538e-04, 2.9126337943686023e-04, 2.9290431507349412e-04;...
                 3.0237666161925968e-04, 3.1368256201851769e-04, 3.1994589783928174e-04, 2.8645692987381806e-04, 2.8398732666288070e-04, 2.8497724812247721e-04;...
                 2.8379960728898217e-04, 2.9493087090087989e-04, 3.0002134460982302e-04, 2.8105264346112700e-04, 2.7992781293530042e-04, 2.8129677837549578e-04;...
                 2.7982528361834261e-04, 2.8149277716683977e-04, 2.8127221745218473e-04, 2.7891422816248457e-04, 2.7764475348614193e-04, 2.7843179501481135e-04;...
                 2.7827560336940959e-04, 2.7742838004891018e-04, 2.7669042671144641e-04, 2.7814103225362878e-04, 2.7725299851353213e-04, 2.7786554070087516e-04;...
                 2.7806430531225439e-04, 2.7728331287889686e-04, 2.7628827947824298e-04, 2.7810727675272061e-04, 2.7725299167447948e-04, 2.7786218181391416e-04;...
                 2.7805307078012390e-04, 2.7725161931786750e-04, 2.7610102877112562e-04, 2.7810727582284800e-04, 2.7725298543407974e-04, 2.7786217396004730e-04;...
                 2.7805083102920585e-04, 2.7724497804197716e-04, 2.7596055271585676e-04, 2.7810727571782779e-04, 2.7725298875309022e-04, 2.7786217335626058e-04;...
                 2.7804759921983664e-04, 2.7722815636335469e-04, 2.7588665848998962e-04, 2.7810727558999303e-04, 2.7725298422004625e-04, 2.7786217242866122e-04];
%                    
supl_pco2i = [   2.7827560336940959e-04, 2.7827560336940959e-04, 2.7827560336940959e-04, 2.7805059701209667e-04, 2.7827560336940959e-04, 2.7827560336940959e-04;... % somblgmt
                 2.7827560336940959e-04, 2.7827560336940959e-04, 2.7827560336940959e-04, 2.7805059701209667e-04, 2.7827560336940959e-04, 2.7827560336940959e-04];     % sombdc

supl_pco2e = [   2.7804812048034909e-04, 2.7725672368938986e-04, 2.7586569242914378e-04, 2.7810767167472760e-04, 2.7725842442262466e-04, 2.7787880134303187e-04;... % somblgmt
                 2.8027693996158463e-04, 2.8537228000959878e-04, 2.8922398685611708e-04, 2.7905276076750055e-04, 2.7781275843732602e-04, 2.7868659733584441e-04];     % sombdc
                        
pco2_anomaly=dustdep_pco2e-dustdep_pco2i;
supl_pco2_anomaly=supl_pco2e-supl_pco2i;
anom_correct=repmat(pco2_anomaly(6,:),[10,1]);
pco2_anomaly=pco2_anomaly-anom_correct;
supl_pco2_anomaly=supl_pco2_anomaly-anom_correct(1:2,:);

%%
%C = linspecer(6,'sequential'); % Colormap for the lines
C=parula(6);

figure
plot(dustdep_frac(:,1),dustdep_pco2e(:,1).*1e6,'o-','color',C(1,:),'MarkerEdgeColor',C(1,:),'MarkerFaceColor',C(1,:),'LineWidth',4)
hold on
plot(dustdep_frac(:,2),dustdep_pco2e(:,2).*1e6,'o-','color',C(2,:),'MarkerEdgeColor',C(2,:),'MarkerFaceColor',C(2,:),'LineWidth',4)
plot(dustdep_frac(:,3),dustdep_pco2e(:,3).*1e6,'o-','color',C(3,:),'MarkerEdgeColor',C(3,:),'MarkerFaceColor',C(3,:),'LineWidth',4)
plot(dustdep_frac(:,4),dustdep_pco2e(:,4).*1e6,'o-','color',C(4,:),'MarkerEdgeColor',C(4,:),'MarkerFaceColor',C(4,:),'LineWidth',4)
plot(dustdep_frac(:,5),dustdep_pco2e(:,5).*1e6,'o-','color',C(5,:),'MarkerEdgeColor',C(5,:),'MarkerFaceColor',C(5,:),'LineWidth',4)
plot(dustdep_frac(:,6),dustdep_pco2e(:,6).*1e6,'o-','color',C(6,:),'MarkerEdgeColor',C(6,:),'MarkerFaceColor',C(6,:),'LineWidth',4)
% Add the spatially varying experiments
scatter(supl_frac(:,1),supl_pco2e(:,1).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(1,:),'LineWidth',1.5)
scatter(supl_frac(:,2),supl_pco2e(:,2).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(2,:),'LineWidth',1.5)
scatter(supl_frac(:,3),supl_pco2e(:,3).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(3,:),'LineWidth',1.5)
scatter(supl_frac(:,4),supl_pco2e(:,4).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(4,:),'LineWidth',1.5)
scatter(supl_frac(:,6),supl_pco2e(:,5).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(5,:),'LineWidth',1.5)
scatter(supl_frac(:,6),supl_pco2e(:,6).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(6,:),'LineWidth',1.5)
legend('1nm lig, a2, dust','1nm lig, a4, dust','1nm lig, a6, dust','1nm lig, a2, all','1nm lig, a4, all','1nm lig, a6, all')
set(gca,'FontSize',16)
xlabel('Ratio of total dust deposition relative to preindustrial','FontSize',16)
ylabel('Atmospheric pCO_2 [uatm]','FontSize',16)
orient landscape
print -dpsc2 ~/Dropbox_Work/PIRES/pires_dust_only_june6.ps

figure
plot(dustdep_frac(:,1),pco2_anomaly(:,1).*1e6,'o-','color',C(1,:),'MarkerEdgeColor',C(1,:),'MarkerFaceColor',C(1,:),'LineWidth',4)
hold on
plot(dustdep_frac(:,2),pco2_anomaly(:,2).*1e6,'o-','color',C(2,:),'MarkerEdgeColor',C(2,:),'MarkerFaceColor',C(2,:),'LineWidth',4)
plot(dustdep_frac(:,3),pco2_anomaly(:,3).*1e6,'o-','color',C(3,:),'MarkerEdgeColor',C(3,:),'MarkerFaceColor',C(3,:),'LineWidth',4)
plot(dustdep_frac(:,4),pco2_anomaly(:,4).*1e6,'o-','color',C(4,:),'MarkerEdgeColor',C(4,:),'MarkerFaceColor',C(4,:),'LineWidth',4)
plot(dustdep_frac(:,5),pco2_anomaly(:,5).*1e6,'o-','color',C(5,:),'MarkerEdgeColor',C(5,:),'MarkerFaceColor',C(5,:),'LineWidth',4)
plot(dustdep_frac(:,6),pco2_anomaly(:,6).*1e6,'o-','color',C(6,:),'MarkerEdgeColor',C(6,:),'MarkerFaceColor',C(6,:),'LineWidth',4)
% Add the spatially varying experiments
scatter(supl_frac(:,1),supl_pco2_anomaly(:,1).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(1,:),'LineWidth',1.5)
scatter(supl_frac(:,2),supl_pco2_anomaly(:,2).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(2,:),'LineWidth',1.5)
scatter(supl_frac(:,3),supl_pco2_anomaly(:,3).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(3,:),'LineWidth',1.5)
scatter(supl_frac(:,4),supl_pco2_anomaly(:,4).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(4,:),'LineWidth',1.5)
scatter(supl_frac(:,5),supl_pco2_anomaly(:,5).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(5,:),'LineWidth',1.5)
scatter(supl_frac(:,6),supl_pco2_anomaly(:,6).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(6,:),'LineWidth',1.5)
legend('1nm lig, a2, dust','1nm lig, a4, dust','1nm lig, a6, dust','1nm lig, a2, all','1nm lig, a4, all','1nm lig, a6, all')
set(gca,'FontSize',16)
xlabel('Ratio of total dust deposition relative to preindustrial','FontSize',16)
ylabel('Atmospheric pCO_2 anomaly [uatm]','FontSize',16)
orient landscape
print -dpsc2 ~/Dropbox_Work/PIRES/pires_dust_only_anomaly_june6.ps

figure
plot(dustdep_mass(:,1),pco2_anomaly(:,1).*1e6,'o-','color',C(1,:),'MarkerEdgeColor',C(1,:),'MarkerFaceColor',C(1,:),'LineWidth',4)
hold on
plot(dustdep_mass(:,2),pco2_anomaly(:,2).*1e6,'o-','color',C(2,:),'MarkerEdgeColor',C(2,:),'MarkerFaceColor',C(2,:),'LineWidth',4)
plot(dustdep_mass(:,3),pco2_anomaly(:,3).*1e6,'o-','color',C(3,:),'MarkerEdgeColor',C(3,:),'MarkerFaceColor',C(3,:),'LineWidth',4)
plot(dustdep_mass(:,4),pco2_anomaly(:,4).*1e6,'o-','color',C(4,:),'MarkerEdgeColor',C(4,:),'MarkerFaceColor',C(4,:),'LineWidth',4)
plot(dustdep_mass(:,5),pco2_anomaly(:,5).*1e6,'o-','color',C(5,:),'MarkerEdgeColor',C(5,:),'MarkerFaceColor',C(5,:),'LineWidth',4)
plot(dustdep_mass(:,6),pco2_anomaly(:,6).*1e6,'o-','color',C(6,:),'MarkerEdgeColor',C(6,:),'MarkerFaceColor',C(6,:),'LineWidth',4)
% Add the spatially varying experiments
scatter(supl_mass(:,1),supl_pco2_anomaly(:,1).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(1,:),'LineWidth',1.5)
scatter(supl_mass(:,2),supl_pco2_anomaly(:,2).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(2,:),'LineWidth',1.5)
scatter(supl_mass(:,3),supl_pco2_anomaly(:,3).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(3,:),'LineWidth',1.5)
scatter(supl_mass(:,4),supl_pco2_anomaly(:,4).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(4,:),'LineWidth',1.5)
scatter(supl_mass(:,5),supl_pco2_anomaly(:,5).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(5,:),'LineWidth',1.5)
scatter(supl_mass(:,6),supl_pco2_anomaly(:,6).*1e6,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(6,:),'LineWidth',1.5)
legend('1nm lig, a2, dust','1nm lig, a4, dust','1nm lig, a6, dust','1nm lig, a2, all','1nm lig, a4, all','1nm lig, a6, all')
set(gca,'FontSize',16)
xlabel('Iron/Dust deposition [Gmol/yr]','FontSize',16)
ylabel('Atmospheric pCO_2 anomaly [uatm]','FontSize',16)
orient landscape
print -dpsc2 ~/Dropbox_Work/PIRES/pires_dust_only_mass_anomaly_june6.ps

%% Ligand and Free iron calculation
ligand_stab=1e8;
freefemax=0.3e-6;
[fein,ligand_tot]=meshgrid([0:1e-7:4.5e-6],[0:1e-7:4e-6]);
lig=(-ligand_stab.*fein + ligand_stab.*ligand_tot-1 + ...
    ((ligand_stab.*fein-ligand_stab.*ligand_tot+1).^2 + ...
    4.*ligand_stab.*ligand_tot).^0.5)/(2.*ligand_stab);
FeL = ligand_tot-lig;
freefe=max(fein-FeL,0);

% adsafe method
thx=freefe;
thy=ones(size(freefe)).*freefemax;
freefe_min_adsafe=( 1 - tanh((thx-thy)./theps) ).* thx./2 +( 1 + tanh((thx-thy)./theps) ).* thy./2;
% not ad friendly, but gives same results
freefe_min=min(freefe,freefemax);

FeT=FeL+freefe_min_adsafe;

figure
contourf(unique(ligand_tot),unique(fein),freefe',[0:1e-7:4.5e-6])
colormap(parula(44));caxis([0 4.5e-6]);colorbar
set(gca,'FontSize',14)
xlabel('Ligand Concentration [mol/m3]','FontSize',14)
ylabel('Iron Concentration [mol/m3]','FontSize',14)
title('Unaltered Free Fe Concentration [mol/m3]','FontSize',14)
orient landscape
print -dpsc2 ~/Desktop/unaltered_free_fe.ps

figure
contourf(unique(ligand_tot),unique(fein),freefe_min',[0:1e-7:4.5e-6])
colormap(parula(44));caxis([0 4.5e-6]);colorbar
set(gca,'FontSize',14)
xlabel('Ligand Concentration [mol/m3]','FontSize',14)
ylabel('Iron Concentration [mol/m3]','FontSize',14)
title('MIN Free Fe Concentration [mol/m3]','FontSize',14)
orient landscape
print -dpsc2 ~/Desktop/min_free_fe.ps

figure
contourf(unique(ligand_tot),unique(fein),FeL',[0:1e-7:4.5e-6])
colormap(parula(44));caxis([0 4.5e-6]);colorbar
set(gca,'FontSize',14)
xlabel('Ligand Concentration [mol/m3]','FontSize',14)
ylabel('Iron Concentration [mol/m3]','FontSize',14)
title('Iron-Ligand Complex Concentration [mol/m3]','FontSize',14)
orient landscape
print -dpsc2 ~/Desktop/conc_iron_ligand_complex.ps

figure
contourf(unique(ligand_tot),unique(fein),FeT',[0:1e-7:4.5e-6])
colormap(parula(44));caxis([0 4.5e-6]);colorbar
set(gca,'FontSize',14)
xlabel('Ligand Concentration [mol/m3]','FontSize',14)
ylabel('Iron Concentration [mol/m3]','FontSize',14)
title('Total Iron Concentration [mol/m3]','FontSize',14)
orient landscape
print -dpsc2 ~/Desktop/conc_total_iron.ps
