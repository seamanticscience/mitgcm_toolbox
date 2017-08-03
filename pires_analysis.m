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
% dustdep_pco2i = [2.7860513774131675,2.8395671452335062,2.7827560336940959,2.7804760178173178;...
%                  2.7860513774131675,2.8395671452335062,2.7827560336940959,2.7804760178173178;...
%                  2.7860513774131675,2.8395671452335062,2.7827560336940959,2.7804760178173178;...
%                  2.7860513774131675,2.8395671452335062,2.7827560336940959,2.7804760178173178;...
%                  2.7827560336940959,2.7827560336940959,2.7827560336940959,2.7827560336940959;...
%                  2.7860513774131675,2.8395671452335062,2.7827560336940959,2.7804760178173178;...
%                  2.7860513774131675,2.8395671452335062,2.7827560336940959,2.7804760178173178;...
%                  2.7860513774131675,2.8395671452335062,2.7827560336940959,2.7804760178173178;...
%                  2.7860513774131675,2.8395671452335062,2.7827560336940959,2.7804760178173178];
% 
% dustdep_pco2e = [2.7865692521628744,2.8737397488010194,3.1602158553843934,3.0391382206044007;...
%                  2.7865070785807459,2.8680560330299849,3.0237666161925968,2.7803668025404266;...
%                  2.7864052600680688,2.8596643595837142,2.8379960728898217,2.7803668025404266;...
%                  2.7862407344460217,2.8505872320980783,2.7982528361834261,2.7803668025404266;...
%                  2.7860513774131675,2.8395671452335062,2.7827560336940959,2.7804760178173178;...
%                  2.7851454017573232,2.8243927127074870,2.7806430531225439,2.7803668025404266;...
%                  2.7835124001088593,2.8072225977773037,2.7805307078012390,2.7803668025404266;...
%                  2.7834765266869770,2.8053255984440429,2.7805083102920585,2.7803668025404266;...
%                  2.7834194770564430,2.8015863208790971,2.7804759921983664,2.7803668025404266];
% 
% %                    
% supl_pco2i = [      2.7860513774131675,	2.8395671452335062,	2.7827560336940959,	2.7804760178173178;... % somblgmt
%                     2.7860513774131675, 2.8395671452335062,	2.7827560336940959,	2.7804760178173178];     % sombdc
% 
% supl_pco2e = [      2.7834416771106371, 2.7804990487806218,	2.7804812048034909,	2.7803668025404266;... % somblgmt
%                     2.7864320877865286, 2.8478539121352595,	2.8027693996158463,	2.7803668025404266];     % sombdc

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

%                        1nm lig          1nm lig            1nm lig               1nm lig            1nm lig             1nm lig  
%                        a2 dust          a4 dust            a6 dust               a2 all             a4 all               a6 all
dustdep_pco2i = [278.27560336940959, 278.27560336940959, 278.27560336940959, 278.05059701209667, 278.27560336940959, 278.27560336940959;...
                 278.27560336940959, 278.27560336940959, 278.27560336940959, 278.05059701209667, 278.27560336940959, 278.27560336940959;...
                 278.27560336940959, 278.27560336940959, 278.27560336940959, 278.05059701209667, 278.27560336940959, 278.27560336940959;...
                 278.27560336940959, 278.27560336940959, 278.27560336940959, 278.05059701209667, 278.27560336940959, 278.27560336940959;...
                 278.27560336940959, 278.27560336940959, 278.27560336940959, 278.05059701209667, 278.27560336940959, 278.27560336940959;...
                 278.27560336940959, 278.27560336940959, 278.27560336940959, 278.05059701209667, 278.27560336940959, 278.27560336940959;...
                 278.27560336940959, 278.27560336940959, 278.27560336940959, 278.05059701209667, 278.27560336940959, 278.27560336940959;...
                 278.27560336940959, 278.27560336940959, 278.27560336940959, 278.05059701209667, 278.27560336940959, 278.27560336940959;...
                 278.27560336940959, 278.27560336940959, 278.27560336940959, 278.05059701209667, 278.27560336940959, 278.27560336940959;...
                 278.27560336940959, 278.27560336940959, 278.27560336940959, 278.05059701209667, 278.27560336940959, 278.27560336940959];

dustdep_pco2e = [317.35176639198438, 340.56542087654746, 346.35531093347229, 298.60463184297470, 294.80437170197233, 294.07540226501739;...
                 316.02158553843934, 336.34911785843758, 343.76981536857480, 297.03380236476538, 291.26337943686023, 292.90431507349412;...
                 302.37666161925968, 313.68256201851769, 319.94589783928174, 286.45692987381806, 283.98732666288070, 284.97724812247721;...
                 283.79960728898217, 294.93087090087989, 300.02134460982302, 281.05264346112700, 279.92781293530042, 281.29677837549578;...
                 279.82528361834261, 281.49277716683977, 281.27221745218473, 278.91422816248457, 277.64475348614193, 278.43179501481135;...
                 278.27560336940959, 277.42838004891018, 276.69042671144641, 278.14103225362878, 277.25299851353213, 277.86554070087516;...
                 278.06430531225439, 277.28331287889686, 276.28827947824298, 278.10727675272061, 277.25299167447948, 277.86218181391416;...
                 278.05307078012390, 277.25161931786750, 276.10102877112562, 278.10727582284800, 277.25298543407974, 277.86217396004730;...
                 278.05083102920585, 277.24497804197716, 275.96055271585676, 278.10727571782779, 277.25298875309022, 277.86217335626058;...
                 278.04759921983664, 277.22815636335469, 275.88665848998962, 278.10727558999303, 277.25298422004625, 277.86217242866122];
%                    
supl_pco2i = [   278.27560336940959, 278.27560336940959, 278.27560336940959, 278.05059701209667, 278.27560336940959, 278.27560336940959;... % somblgmt
                 278.27560336940959, 278.27560336940959, 278.27560336940959, 278.05059701209667, 278.27560336940959, 278.27560336940959];     % sombdc

supl_pco2e = [   278.04812048034909, 277.25672368938986, 275.86569242914378, 278.10767167472760, 277.25842442262466, 277.87880134303187;... % somblgmt
                 280.27693996158463, 285.37228000959878, 289.22398685611708, 279.05276076750055, 277.81275843732602, 278.68659733584441];     % sombdc

%                       0.25nm              1nm                 4nm
%                       a2 dust           a2 dust            a2 dust
ligand_pco2i = [ 280.34506024059386, 278.27560336940959, 278.04760178173178;...
                 280.34506024059386, 278.27560336940959, 278.04760178173178;...
                 280.34506024059386, 278.27560336940959, 278.04760178173178;...
                 280.34506024059386, 278.27560336940959, 278.04760178173178;...
                 280.34506024059386, 278.27560336940959, 278.04760178173178;...
                 280.34506024059386, 278.27560336940959, 278.27560336940959;...
                 280.34506024059386, 278.27560336940959, 278.04760178173178;...
                 280.34506024059386, 278.27560336940959, 278.04760178173178;...
                 280.34506024059386, 278.27560336940959, 278.04760178173178;...
                 280.34506024059386, 278.27560336940959, 278.04760178173178 ];
            
ligand_pco2e = [ 306.13403846319021, 317.35176639198438, 373.14287896779197;...
                 305.60143653266290, 316.02158553843934, 345.40377281355527;...
                 300.59482096263197, 302.37666161925968, 278.07519692254578;...
                 293.62160187391322, 283.79960728898217, 278.02678624345468;...
                 287.69253188154503, 279.82528361834261, 278.03668025404266;...
                 280.32130277419221, 278.27560336940959, 278.04760178173178;...
                 271.43113237032264, 278.06430531225439, 278.03668025404266;...
                 263.94755001404315, 278.05307078012390, 278.03668025404266;...
                 261.93326164913096, 278.05083102920585, 278.03668025404266;...
                 255.90376134088224, 278.04759921983664, 278.03668025404266 ];

             
supl_lig_pco2i=[ 280.34506024059386, 278.27560336940959, 278.04760178173178;...
                 280.34506024059386, 278.27560336940959, 278.04760178173178 ];
             
supl_lig_pco2e=[ 249.48588242902690, 278.04812048034909, 278.03668025404266;...
                 286.97632288776729, 280.27693996158463, 278.03668025404266 ];

%                 1nm lig     1nm lig   1nm lig     1nm lig   1nm lig   1nm lig
%                 a2 dust     a4 dust   a6 dust     a2 all    a4 all    a6 all
dustdep_pstar =[  0.2628,     0.3417,   0.3939,     0.3219,   0.4754,   0.5344;...
                  0.2668,     0.3521,   0.4004,     0.3271,   0.4852,   0.5361;...
                  0.3090,     0.4152,   0.4613,     0.3633,   0.5038,   0.5484;...
                  0.3716,     0.4738,   0.5211,     0.3825,   0.5125,   0.5532;...
                  0.3858,     0.5094,   0.5502,     0.3901,   0.5167,   0.5552;...
                  0.3910,     0.5171,   0.5540,     0.3928,   0.5174,   0.5553;...
                  0.3918,     0.5173,   0.5545,     0.3929,   0.5174,   0.5553;...
                  0.3918,     0.5174,   0.5548,     0.3929,   0.5174,   0.5553;...
                  0.3918,     0.5174,   0.5550,     0.3929,   0.5174,   0.5553;...
                  0.3919,     0.5174,   0.5552,     0.3929,   0.5174,   0.5553 ];


supl_pstar  =  [  0.3919,     0.5174,   0.5554,     0.3929,   0.5174,   0.5553;...
                  0.3846,     0.5007,   0.5384,     0.3899,   0.5163    0.5539 ];       

ligand_pstar  =[  0.1937,     0.2628,   0.1307;...
                  0.1954,     0.2668,   0.1977;...
                  0.2112,     0.3090,   0.3917;...
                  0.2340,     0.3716,   0.3919;...
                  0.2540,     0.3858,   0.3920;...
                  0.2797,     0.3910,   0.3919;...
                  0.3111,     0.3918,   0.3920;...
                  0.3377,     0.3918,   0.3920;...
                  0.3447,     0.3918,   0.3920;...
                  0.3654,     0.3919,   0.3920 ];


supl_lig_pstar=[ 0.3849,      0.3919,   0.3920;...
                 0.2572,      0.3846,   0.3920 ];

             
dustdep_bp   = [ 15.9918,     17.5033,   18.4128,   20.0698,   26.4333,   27.7891;...
                 16.3004,     18.2707,   18.8402,   20.4876,   27.5557,   28.1706;...
                 19.5545,     23.3694,   23.4194,   23.3228,   29.6967,   30.5831;...
                 24.0752,     27.3363,   27.7337,   24.8625,   30.7625,   31.3711;...
                 25.3403,     30.3904,   31.1544,   25.6817,   31.6829,   32.2542;...
                 25.9334,     31.8631,   32.5099,   26.0484,   31.9260,   32.5533;...
                 26.0642,     31.9238,   32.5467,   26.0776,   31.9260,   32.5579;...
                 26.0646,     31.9246,   32.5507,   26.0776,   31.9260,   32.5579;...
                 26.0647,     31.9248,   32.5537,   26.0776,   31.9260,   32.5579;...
                 26.0648,     31.9254,   32.5554,   26.0776,   31.9260,   32.5579];
    

supl_bp      = [ 26.0642,     31.8908,   32.5043,   26.0772,   31.9179,   32.5365;...
                 25.2269,     29.4635,   29.4186,   25.6356,   31.5962,   32.1996];



ligand_bp    = [ 13.2159,     15.9918,   6.2354 ;...
                 13.3166,     16.3004,   10.7228;...
                 14.2860,     19.5545,   26.0185;...
                 15.7631,     24.0752,   26.0609;...
                 17.1680,     25.3403,   26.0635;...
                 19.0175,     25.9334,   26.0650;...
                 21.3610,     26.0642,   26.0635;...
                 23.4453,     26.0646,   26.0635;...
                 23.9978,     26.0647,   26.0635;...
                 25.4240,     26.0648,   26.0635];
                 
supl_lig_bp  = [ 25.4282,     26.0642,   26.0635;...
                 17.3818,     25.2269,   26.0635];
    
% For 1 nm ligand experiments
pco2_anomaly=dustdep_pco2e-dustdep_pco2i;
supl_pco2_anomaly=supl_pco2e-supl_pco2i;
anom_correct=repmat(pco2_anomaly(6,:),[10,1]);
pco2_anomaly=pco2_anomaly-anom_correct;
supl_pco2_anomaly=supl_pco2_anomaly-anom_correct(1:2,:);

% For varying ligand experiments
lig_pco2_anomaly=ligand_pco2e-ligand_pco2i;
supl_lig_pco2_anomaly=supl_lig_pco2e-supl_lig_pco2i;
anom_correct=repmat(lig_pco2_anomaly(6,:),[10,1]);
lig_pco2_anomaly=lig_pco2_anomaly-anom_correct;
supl_lig_pco2_anomaly=supl_lig_pco2_anomaly-anom_correct(1:2,:);

% Add up the outgassing or uptake during initial control run from the
% change in DIC concentration to find change in pCO2 if it had not been
% held constant.
lig_pco2=ligand_pco2e;
lig_pco2(:,1)=lig_pco2(:,1)+181;
lig_pco2(:,3)=lig_pco2(:,3)-0.228 ;
supl_lig_pco2=supl_lig_pco2e;
supl_lig_pco2(:,1)=supl_lig_pco2(:,1)+181;
supl_lig_pco2(:,3)=supl_lig_pco2(:,3)-0.228;

%%
%C = linspecer(6,'sequential'); % Colormap for the lines
C=parula(6);

% figure
% plot(dustdep_frac(:,1),dustdep_pco2e(:,1),'o-','color',C(1,:),'MarkerEdgeColor',C(1,:),'MarkerFaceColor',C(1,:),'LineWidth',4)
% hold on
% plot(dustdep_frac(:,2),dustdep_pco2e(:,2),'o-','color',C(2,:),'MarkerEdgeColor',C(2,:),'MarkerFaceColor',C(2,:),'LineWidth',4)
% plot(dustdep_frac(:,3),dustdep_pco2e(:,3),'o-','color',C(3,:),'MarkerEdgeColor',C(3,:),'MarkerFaceColor',C(3,:),'LineWidth',4)
% plot(dustdep_frac(:,4),dustdep_pco2e(:,4),'o-','color',C(4,:),'MarkerEdgeColor',C(4,:),'MarkerFaceColor',C(4,:),'LineWidth',4)
% plot(dustdep_frac(:,5),dustdep_pco2e(:,5),'o-','color',C(5,:),'MarkerEdgeColor',C(5,:),'MarkerFaceColor',C(5,:),'LineWidth',4)
% plot(dustdep_frac(:,6),dustdep_pco2e(:,6),'o-','color',C(6,:),'MarkerEdgeColor',C(6,:),'MarkerFaceColor',C(6,:),'LineWidth',4)
% % Add the spatially varying experiments
% scatter(supl_frac(:,1),supl_pco2e(:,1),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(1,:),'LineWidth',1.5)
% scatter(supl_frac(:,2),supl_pco2e(:,2),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(2,:),'LineWidth',1.5)
% scatter(supl_frac(:,3),supl_pco2e(:,3),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(3,:),'LineWidth',1.5)
% scatter(supl_frac(:,4),supl_pco2e(:,4),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(4,:),'LineWidth',1.5)
% scatter(supl_frac(:,6),supl_pco2e(:,5),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(5,:),'LineWidth',1.5)
% scatter(supl_frac(:,6),supl_pco2e(:,6),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(6,:),'LineWidth',1.5)
% legend('1nm lig, a2, dust','1nm lig, a4, dust','1nm lig, a6, dust','1nm lig, a2, all','1nm lig, a4, all','1nm lig, a6, all')
% set(gca,'FontSize',16)
% xlabel('Ratio of total dust deposition relative to preindustrial','FontSize',16)
% ylabel('Atmospheric pCO_2 [uatm]','FontSize',16)
% orient landscape
% print('-dpsc2',['~/Dropbox_Work/PIRES/pires_dust_perturb_',date,'.ps']);

figure
plot(dustdep_frac(:,1),pco2_anomaly(:,1),'o-','color',C(1,:),'MarkerEdgeColor',C(1,:),'MarkerFaceColor',C(1,:),'LineWidth',4)
hold on
plot(dustdep_frac(:,2),pco2_anomaly(:,2),'o-','color',C(2,:),'MarkerEdgeColor',C(2,:),'MarkerFaceColor',C(2,:),'LineWidth',4)
plot(dustdep_frac(:,3),pco2_anomaly(:,3),'o-','color',C(3,:),'MarkerEdgeColor',C(3,:),'MarkerFaceColor',C(3,:),'LineWidth',4)
plot(dustdep_frac(:,4),pco2_anomaly(:,4),'o-','color',C(4,:),'MarkerEdgeColor',C(4,:),'MarkerFaceColor',C(4,:),'LineWidth',4)
plot(dustdep_frac(:,5),pco2_anomaly(:,5),'o-','color',C(5,:),'MarkerEdgeColor',C(5,:),'MarkerFaceColor',C(5,:),'LineWidth',4)
plot(dustdep_frac(:,6),pco2_anomaly(:,6),'o-','color',C(6,:),'MarkerEdgeColor',C(6,:),'MarkerFaceColor',C(6,:),'LineWidth',4)
% Add the spatially varying experiments
scatter(supl_frac(:,1),supl_pco2_anomaly(:,1),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(1,:),'LineWidth',1.5)
scatter(supl_frac(:,2),supl_pco2_anomaly(:,2),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(2,:),'LineWidth',1.5)
scatter(supl_frac(:,3),supl_pco2_anomaly(:,3),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(3,:),'LineWidth',1.5)
scatter(supl_frac(:,4),supl_pco2_anomaly(:,4),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(4,:),'LineWidth',1.5)
scatter(supl_frac(:,5),supl_pco2_anomaly(:,5),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(5,:),'LineWidth',1.5)
scatter(supl_frac(:,6),supl_pco2_anomaly(:,6),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(6,:),'LineWidth',1.5)
legend('1nm lig, a2, dust','1nm lig, a4, dust','1nm lig, a6, dust','1nm lig, a2, all','1nm lig, a4, all','1nm lig, a6, all')
set(gca,'FontSize',16,'XLim',[0 5],'XTick',1:1:5)
xlabel('Ratio of total dust deposition relative to preindustrial','FontSize',16)
ylabel('Atmospheric pCO_2 anomaly [uatm]','FontSize',16)
orient landscape
print('-dpsc2',['~/Dropbox_Work/PIRES/pires_dust_perturb_anomaly_',date,'.ps']);

% figure
% plot(dustdep_mass(:,1),pco2_anomaly(:,1),'o-','color',C(1,:),'MarkerEdgeColor',C(1,:),'MarkerFaceColor',C(1,:),'LineWidth',4)
% hold on
% plot(dustdep_mass(:,2),pco2_anomaly(:,2),'o-','color',C(2,:),'MarkerEdgeColor',C(2,:),'MarkerFaceColor',C(2,:),'LineWidth',4)
% plot(dustdep_mass(:,3),pco2_anomaly(:,3),'o-','color',C(3,:),'MarkerEdgeColor',C(3,:),'MarkerFaceColor',C(3,:),'LineWidth',4)
% plot(dustdep_mass(:,4),pco2_anomaly(:,4),'o-','color',C(4,:),'MarkerEdgeColor',C(4,:),'MarkerFaceColor',C(4,:),'LineWidth',4)
% plot(dustdep_mass(:,5),pco2_anomaly(:,5),'o-','color',C(5,:),'MarkerEdgeColor',C(5,:),'MarkerFaceColor',C(5,:),'LineWidth',4)
% plot(dustdep_mass(:,6),pco2_anomaly(:,6),'o-','color',C(6,:),'MarkerEdgeColor',C(6,:),'MarkerFaceColor',C(6,:),'LineWidth',4)
% % Add the spatially varying experiments
% scatter(supl_mass(:,1),supl_pco2_anomaly(:,1),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(1,:),'LineWidth',1.5)
% scatter(supl_mass(:,2),supl_pco2_anomaly(:,2),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(2,:),'LineWidth',1.5)
% scatter(supl_mass(:,3),supl_pco2_anomaly(:,3),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(3,:),'LineWidth',1.5)
% scatter(supl_mass(:,4),supl_pco2_anomaly(:,4),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(4,:),'LineWidth',1.5)
% scatter(supl_mass(:,5),supl_pco2_anomaly(:,5),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(5,:),'LineWidth',1.5)
% scatter(supl_mass(:,6),supl_pco2_anomaly(:,6),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(6,:),'LineWidth',1.5)
% legend('1nm lig, a2, dust','1nm lig, a4, dust','1nm lig, a6, dust','1nm lig, a2, all','1nm lig, a4, all','1nm lig, a6, all')
% set(gca,'FontSize',16)
% xlabel('Iron/Dust deposition [Gmol/yr]','FontSize',16)
% ylabel('Atmospheric pCO_2 anomaly [uatm]','FontSize',16)
% orient landscape
% print('-dpsc2',['~/Dropbox_Work/PIRES/pires_dust_perturb_mass_anomaly_',date,'.ps']);

figure
plot(dustdep_frac(:,1),lig_pco2_anomaly(:,1),'o-','color',C(1,:),'MarkerEdgeColor',C(1,:),'MarkerFaceColor',C(1,:),'LineWidth',4)
hold on
plot(dustdep_frac(:,2),lig_pco2_anomaly(:,2),'o-','color',C(4,:),'MarkerEdgeColor',C(4,:),'MarkerFaceColor',C(4,:),'LineWidth',4)
plot(dustdep_frac(:,3),lig_pco2_anomaly(:,3),'o-','color',C(5,:),'MarkerEdgeColor',C(5,:),'MarkerFaceColor',C(5,:),'LineWidth',4)
%plot(dustdep_frac(:,4),lig_pco2_anomaly(:,4),'o-','color',C(4,:),'MarkerEdgeColor',C(4,:),'MarkerFaceColor',C(4,:),'LineWidth',4)
%plot(dustdep_frac(:,5),lig_pco2_anomaly(:,5),'o-','color',C(5,:),'MarkerEdgeColor',C(5,:),'MarkerFaceColor',C(5,:),'LineWidth',4)
%plot(dustdep_frac(:,6),lig_pco2_anomaly(:,6),'o-','color',C(6,:),'MarkerEdgeColor',C(6,:),'MarkerFaceColor',C(6,:),'LineWidth',4)
% Add the spatially varying experiments
scatter(supl_frac(:,1),supl_lig_pco2_anomaly(:,1),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(1,:),'LineWidth',1.5)
scatter(supl_frac(:,2),supl_lig_pco2_anomaly(:,2),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(4,:),'LineWidth',1.5)
scatter(supl_frac(:,3),supl_lig_pco2_anomaly(:,3),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(5,:),'LineWidth',1.5)
%scatter(supl_frac(:,4),supl_lig_pco2_anomaly(:,4),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(4,:),'LineWidth',1.5)
%scatter(supl_frac(:,5),supl_lig_pco2_anomaly(:,5),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(5,:),'LineWidth',1.5)
%scatter(supl_frac(:,6),supl_lig_pco2_anomaly(:,6),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(6,:),'LineWidth',1.5)
legend('0.25nm lig, a2, dust','1nm lig, a2, dust','4nm lig, a2, dust') %,'1nm lig, a2, all','1nm lig, a4, all','1nm lig, a6, all')
set(gca,'FontSize',16,'XLim',[0 5],'XTick',1:1:5)
xlabel('Ratio of total dust deposition relative to preindustrial','FontSize',16)
ylabel('Atmospheric pCO_2 anomaly [uatm]','FontSize',16)
orient landscape
print('-dpsc2',['~/Dropbox_Work/PIRES/pires_ligand_dust_perturb_anomaly_',date,'.ps']);

figure
plot(dustdep_frac(:,1),lig_pco2(:,1),'o-','color',C(1,:),'MarkerEdgeColor',C(1,:),'MarkerFaceColor',C(1,:),'LineWidth',4)
hold on
plot(dustdep_frac(:,2),lig_pco2(:,2),'o-','color',C(4,:),'MarkerEdgeColor',C(4,:),'MarkerFaceColor',C(4,:),'LineWidth',4)
plot(dustdep_frac(:,3),lig_pco2(:,3),'o-','color',C(5,:),'MarkerEdgeColor',C(5,:),'MarkerFaceColor',C(5,:),'LineWidth',4)
% plot(dustdep_frac(:,4),lig_pco2(:,4),'o-','color',C(4,:),'MarkerEdgeColor',C(4,:),'MarkerFaceColor',C(4,:),'LineWidth',4)
% plot(dustdep_frac(:,5),lig_pco2(:,5),'o-','color',C(5,:),'MarkerEdgeColor',C(5,:),'MarkerFaceColor',C(5,:),'LineWidth',4)
% plot(dustdep_frac(:,6),lig_pco2(:,6),'o-','color',C(6,:),'MarkerEdgeColor',C(6,:),'MarkerFaceColor',C(6,:),'LineWidth',4)
% Add the spatially varying experiments
scatter(supl_frac(:,1),supl_lig_pco2(:,1),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(1,:),'LineWidth',1.5)
scatter(supl_frac(:,2),supl_lig_pco2(:,2),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(2,:),'LineWidth',1.5)
scatter(supl_frac(:,3),supl_lig_pco2(:,3),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(4,:),'LineWidth',1.5)
% scatter(supl_frac(:,4),supl_lig_pco2(:,4),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(4,:),'LineWidth',1.5)
% scatter(supl_frac(:,6),supl_lig_pco2(:,5),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(5,:),'LineWidth',1.5)
% scatter(supl_frac(:,6),supl_lig_pco2(:,6),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(6,:),'LineWidth',1.5)
legend('0.25nm lig, a2, dust','1nm lig, a2, dust','4nm lig, a2, dust') %,'1nm lig, a2, all','1nm lig, a4, all','1nm lig, a6, all')
set(gca,'FontSize',16,'XLim',[0 5],'XTick',1:1:5)
xlabel('Ratio of total dust deposition relative to preindustrial','FontSize',16)
ylabel('Atmospheric pCO_2 [uatm]','FontSize',16)
orient landscape
print('-dpsc2',['~/Dropbox_Work/PIRES/pires_ligand_dust_perturb_',date,'.ps']);

%%
figure
plot(dustdep_frac(:,1),dustdep_pstar(:,1).*1e2,'o-','color',C(1,:),'MarkerEdgeColor',C(1,:),'MarkerFaceColor',C(1,:),'LineWidth',4)
hold on
plot(dustdep_frac(:,2),dustdep_pstar(:,2).*1e2,'o-','color',C(2,:),'MarkerEdgeColor',C(2,:),'MarkerFaceColor',C(2,:),'LineWidth',4)
plot(dustdep_frac(:,3),dustdep_pstar(:,3).*1e2,'o-','color',C(3,:),'MarkerEdgeColor',C(3,:),'MarkerFaceColor',C(3,:),'LineWidth',4)
plot(dustdep_frac(:,4),dustdep_pstar(:,4).*1e2,'o-','color',C(4,:),'MarkerEdgeColor',C(4,:),'MarkerFaceColor',C(4,:),'LineWidth',4)
plot(dustdep_frac(:,5),dustdep_pstar(:,5).*1e2,'o-','color',C(5,:),'MarkerEdgeColor',C(5,:),'MarkerFaceColor',C(5,:),'LineWidth',4)
plot(dustdep_frac(:,6),dustdep_pstar(:,6).*1e2,'o-','color',C(6,:),'MarkerEdgeColor',C(6,:),'MarkerFaceColor',C(6,:),'LineWidth',4)
% Add the spatially varying experiments
scatter(supl_frac(:,1),supl_pstar(:,1).*1e2,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(1,:),'LineWidth',1.5)
scatter(supl_frac(:,2),supl_pstar(:,2).*1e2,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(2,:),'LineWidth',1.5)
scatter(supl_frac(:,3),supl_pstar(:,3).*1e2,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(3,:),'LineWidth',1.5)
scatter(supl_frac(:,4),supl_pstar(:,4).*1e2,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(4,:),'LineWidth',1.5)
scatter(supl_frac(:,5),supl_pstar(:,5).*1e2,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(5,:),'LineWidth',1.5)
scatter(supl_frac(:,6),supl_pstar(:,6).*1e2,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(6,:),'LineWidth',1.5)
legend('1nm lig, a2, dust','1nm lig, a4, dust','1nm lig, a6, dust','1nm lig, a2, all','1nm lig, a4, all','1nm lig, a6, all','Location','SouthEast')
set(gca,'FontSize',16,'XLim',[0 5],'XTick',1:1:5)
xlabel('Ratio of total dust deposition relative to preindustrial','FontSize',16)
ylabel('Pstar metric [%]','FontSize',16)
orient landscape
print('-dpsc2',['~/Dropbox_Work/PIRES/pires_dust_perturb_pstar_',date,'.ps']);

figure
plot(dustdep_frac(:,1),ligand_pstar(:,1).*1e2,'o-','color',C(1,:),'MarkerEdgeColor',C(1,:),'MarkerFaceColor',C(1,:),'LineWidth',4)
hold on
plot(dustdep_frac(:,2),ligand_pstar(:,2).*1e2,'o-','color',C(4,:),'MarkerEdgeColor',C(4,:),'MarkerFaceColor',C(4,:),'LineWidth',4)
plot(dustdep_frac(:,3),ligand_pstar(:,3).*1e2,'o-','color',C(5,:),'MarkerEdgeColor',C(5,:),'MarkerFaceColor',C(5,:),'LineWidth',4)
%plot(dustdep_frac(:,4),ligand_pstar(:,4).*1e2,'o-','color',C(4,:),'MarkerEdgeColor',C(4,:),'MarkerFaceColor',C(4,:),'LineWidth',4)
%plot(dustdep_frac(:,5),ligand_pstar(:,5).*1e2,'o-','color',C(5,:),'MarkerEdgeColor',C(5,:),'MarkerFaceColor',C(5,:),'LineWidth',4)
%plot(dustdep_frac(:,6),ligand_pstar(:,6).*1e2,'o-','color',C(6,:),'MarkerEdgeColor',C(6,:),'MarkerFaceColor',C(6,:),'LineWidth',4)
% Add the spatially varying experiments
scatter(supl_frac(:,1),supl_lig_pstar(:,1).*1e2,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(1,:),'LineWidth',1.5)
scatter(supl_frac(:,2),supl_lig_pstar(:,2).*1e2,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(4,:),'LineWidth',1.5)
scatter(supl_frac(:,3),supl_lig_pstar(:,3).*1e2,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(5,:),'LineWidth',1.5)
%scatter(supl_frac(:,4),supl_lig_pstar(:,4).*1e2,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(4,:),'LineWidth',1.5)
%scatter(supl_frac(:,5),supl_lig_pstar(:,5).*1e2,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(5,:),'LineWidth',1.5)
%scatter(supl_frac(:,6),supl_lig_pstar(:,6).*1e2,500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(6,:),'LineWidth',1.5)
legend('0.25nm lig, a2, dust','1nm lig, a2, dust','4nm lig, a2, dust','Location','SouthEast') %,'1nm lig, a2, all','1nm lig, a4, all','1nm lig, a6, all')
set(gca,'FontSize',16,'XLim',[0 5],'XTick',1:1:5)
xlabel('Ratio of total dust deposition relative to preindustrial','FontSize',16)
ylabel('Pstar metric [%]','FontSize',16)
orient landscape
print('-dpsc2',['~/Dropbox_Work/PIRES/pires_ligand_dust_perturb_pstar_',date,'.ps']);

%%
figure
plot(dustdep_frac(:,1),dustdep_bp(:,1),'o-','color',C(1,:),'MarkerEdgeColor',C(1,:),'MarkerFaceColor',C(1,:),'LineWidth',4)
hold on
plot(dustdep_frac(:,2),dustdep_bp(:,2),'o-','color',C(2,:),'MarkerEdgeColor',C(2,:),'MarkerFaceColor',C(2,:),'LineWidth',4)
plot(dustdep_frac(:,3),dustdep_bp(:,3),'o-','color',C(3,:),'MarkerEdgeColor',C(3,:),'MarkerFaceColor',C(3,:),'LineWidth',4)
plot(dustdep_frac(:,4),dustdep_bp(:,4),'o-','color',C(4,:),'MarkerEdgeColor',C(4,:),'MarkerFaceColor',C(4,:),'LineWidth',4)
plot(dustdep_frac(:,5),dustdep_bp(:,5),'o-','color',C(5,:),'MarkerEdgeColor',C(5,:),'MarkerFaceColor',C(5,:),'LineWidth',4)
plot(dustdep_frac(:,6),dustdep_bp(:,6),'o-','color',C(6,:),'MarkerEdgeColor',C(6,:),'MarkerFaceColor',C(6,:),'LineWidth',4)
% Add the spatially varying experiments
scatter(supl_frac(:,1),supl_bp(:,1),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(1,:),'LineWidth',1.5)
scatter(supl_frac(:,2),supl_bp(:,2),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(2,:),'LineWidth',1.5)
scatter(supl_frac(:,3),supl_bp(:,3),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(3,:),'LineWidth',1.5)
scatter(supl_frac(:,4),supl_bp(:,4),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(4,:),'LineWidth',1.5)
scatter(supl_frac(:,5),supl_bp(:,5),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(5,:),'LineWidth',1.5)
scatter(supl_frac(:,6),supl_bp(:,6),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(6,:),'LineWidth',1.5)
legend('1nm lig, a2, dust','1nm lig, a4, dust','1nm lig, a6, dust','1nm lig, a2, all','1nm lig, a4, all','1nm lig, a6, all','Location','SouthEast')
set(gca,'FontSize',16,'XLim',[0 5],'XTick',1:1:5)
xlabel('Ratio of total dust deposition relative to preindustrial','FontSize',16)
ylabel('Integrated Biological Production [GtC]','FontSize',16)
orient landscape
print('-dpsc2',['~/Dropbox_Work/PIRES/pires_dust_perturb_bp_',date,'.ps']);

figure
plot(dustdep_frac(:,1),ligand_bp(:,1),'o-','color',C(1,:),'MarkerEdgeColor',C(1,:),'MarkerFaceColor',C(1,:),'LineWidth',4)
hold on
plot(dustdep_frac(:,2),ligand_bp(:,2),'o-','color',C(4,:),'MarkerEdgeColor',C(4,:),'MarkerFaceColor',C(4,:),'LineWidth',4)
plot(dustdep_frac(:,3),ligand_bp(:,3),'o-','color',C(5,:),'MarkerEdgeColor',C(5,:),'MarkerFaceColor',C(5,:),'LineWidth',4)
%plot(dustdep_frac(:,4),ligand_bp(:,4),'o-','color',C(4,:),'MarkerEdgeColor',C(4,:),'MarkerFaceColor',C(4,:),'LineWidth',4)
%plot(dustdep_frac(:,5),ligand_bp(:,5),'o-','color',C(5,:),'MarkerEdgeColor',C(5,:),'MarkerFaceColor',C(5,:),'LineWidth',4)
%plot(dustdep_frac(:,6),ligand_bp(:,6),'o-','color',C(6,:),'MarkerEdgeColor',C(6,:),'MarkerFaceColor',C(6,:),'LineWidth',4)
% Add the spatially varying experiments
scatter(supl_frac(:,1),supl_lig_bp(:,1),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(1,:),'LineWidth',1.5)
scatter(supl_frac(:,2),supl_lig_bp(:,2),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(4,:),'LineWidth',1.5)
scatter(supl_frac(:,3),supl_lig_bp(:,3),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(5,:),'LineWidth',1.5)
%scatter(supl_frac(:,4),supl_lig_bp(:,4),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(4,:),'LineWidth',1.5)
%scatter(supl_frac(:,5),supl_lig_bp(:,5),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(5,:),'LineWidth',1.5)
%scatter(supl_frac(:,6),supl_lig_bp(:,6),500,'p','filled','MarkerEdgeColor','k','MarkerFaceColor',C(6,:),'LineWidth',1.5)
legend('0.25nm lig, a2, dust','1nm lig, a2, dust','4nm lig, a2, dust','Location','SouthEast') %,'1nm lig, a2, all','1nm lig, a4, all','1nm lig, a6, all')
set(gca,'FontSize',16,'XLim',[0 5],'XTick',1:1:5)
xlabel('Ratio of total dust deposition relative to preindustrial','FontSize',16)
ylabel('Integrated Biological Production [GtC]','FontSize',16)
orient landscape
print('-dpsc2',['~/Dropbox_Work/PIRES/pires_ligand_dust_perturb_bp_',date,'.ps']);

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
