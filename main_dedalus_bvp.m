clear all;
close all;
% slurm_num={'12516590',
%     '12516677',
%     '12516679',
%     '12519906',
%     '12520491',
%     '12632332',
%     '12632508',
%     '12632582',
%     '12635575',
%     '12639413'};

% group_name='HB_porous_Nu_kx_Ra';
% group_name='HB_porous_Nu_kx_Ra';
% group_name='trevisan_contour';
% group_name='rosenberg_time_dependent';
% group_name='HB_benard_yang';
% group_name='HB_benard_diffusive_shear';
% group_name='hewitt_2_layer_Omega';
% group_name='hewitt_2_layer_Omega';
% group_name='HB_porous_kx';
% group_name='HB_benard_salt_finger_R_rho_T2S';
% group_name='HB_benard_salt_finger_tau';
% group_name='HB_benard_salt_finger_kx'
% group_name='HB_benard_salt_finger_profile_S1_yang';
% group_name = 'HB_benard_salt_finger_profile';
% group_name='HB_benard_diffusive_kx';
%All of these are for porous media
% group_name='HB_porous_mamou_Ra_T_53_SW';
group_name='HB_porous_rosenberg_Ra_T_100_SW';
% group_name='HB_porous_mamou_Ra_T_53_SW';
% group_name='HB_porous_rosenberg_Ra_T_100_SW';
switch group_name
    case 'HB_porous_thermal_BC'
        %These are results comparing only for the thermal convection
        slurm_num={'12639319',
                '12640193',
                '12640192',
                '12639413',
                '12640102'
                }; %'12640104',
    case 'HB_porous_diffusive_BC'
        %These are doing results for double-diffusive convection... 
        slurm_num={'12643033',
                    '12643076',
                    '12643064',
                    '12643209',
                    '12643124',
                    '12643067',
                    '12643068'};
    case 'HB_porous_finger_BC'
        % These for salt-finger regime...
        slurm_num={'12644479',
                   '12646311',
                   '12704259',
                   '12646317',
                   '12646319',
                   '12646320',
                   '12646324'}; %           '12646314',
    case 'HB_porous_DDC'
        slurm_num={'12997946'};
        
    case 'HB_porous_kx'
        slurm_num={'12784649'};
        %slurm_num={'12761570'}%,...%This is continued from 12639319
                   %'12761797'}; %This is continued from the previous case. 
    case 'HB_porous_Nu_kx_Ra'
        slurm_num={'12784649',...
                    '12784650',...
                    '12784652',...
                    '12784653',...
                    '12784654',...
                    '12784655',...
                    '12784656',...
                    '12784658',...
                    '12784660',...
                    '12784661'};
%           slurm_num=slurm_num(1);
    case 'HB_porous_IVP'
        %slurm_num={'12989358','12989675','12989678'};
        slurm_num={'12991007','12991008','12991009'};
    case 'hewitt_2D'
        %%These three case is comparing with the z depedence of the mean
        %%temeperature and rms..
        slurm_num={'12687510',
                    '12687731',
                    '12687732'}; %Hewitt et al. (2012) 2D
    case 'hewitt_2D_Ra'
        %%This result is comparing with the Ra trend of Nu number, rms..
        slurm_num={'12822713'};
    case 'hewitt_3D'
        %%These three case is comparing with the z depedence of the mean
        %%temeperature and rms..
        slurm_num={'12687744',
                  '12687755',
                  '12687756'}; %Hewitt et al. (2014) 3D
%         flag.uvw_hewitt=1;
              %'12640113', %%This case is wierd... not right
    
    case 'hewitt_3D_Ra'
        slurm_num={'12823940'};
    case 'wen_chini_2D'
        slurm_num={'12784615'};
        
    case 'HB_benard_kx'
        slurm_num={'12760884'}; %Ra=10^6, kx=1~20, reproducing figure 7 of Toomre (1977)
    case 'HB_benard_Ra'
        slurm_num={'12761067'}; %kx=1, Ra=10^6~10^10, reproducing figure 4 of Toomre (1977)
    case 'HB_benard_yang'
%         slurm_num={'13100451', %, 
%                    '13100497'};
         slurm_num={'13106633'};
%         slurm_num=slurm_num(1);
    case 'HB_benard_salt_finger_test'
        slurm_num={'13107994',%R_\rho=Ra_T/Ra_S=2
                   '13109064',%%density ratio R_\rho=Ra_T/Ra_S=5
                   '13109065',%R_\rho=10
                   '13109799',%R_\rho=50
                   '13109795',%R_\rho=2, but neumann boundary condition for u and v velocity
                   '13145554',%R_\rho=2, but with tau=0.1. 
                   '13113091',%R_\rho=2, but with tau=1/3 
                   '13146399',%%tau continuation, from tau=0.01,0.02,0.03,0.04,0.05  %%%0.03,1/3,0.03,0.01
                   '13146406',%R_rho_S2T continuation, with R_rho_S2T=2,3,4...50
                   '13146408',%R_rho_S2T=2, tau=0.1, W0=1 initial condition only 1 layer
                   '13146411',%R_rho_S2T=2, tau=0.1, W0=10, 100, initial condition, solution show the three layer structure
                   '13146434',%R_rho_S2T=2, tau=0.1, W0=1000, shows asymmetric solution...
                   '13178474',%R_rho_S2T=2, tau=0.1, asymmetricsolution, continued from the previous one.
                   '13207297',%IVP, 13146408+noise
                   '13207298',%IVP, 13146411+noise
                   '13207299'%IVP, 13146434+noise
                   %'13207282',%IVP, %R_rho_S2T=2, tau=0.1, W0=1 initial condition only 1 layer
                   %'13207283',%IVP, %R_rho_S2T=2, tau=0.1, W0=10 initial condition only 1 layer
                   %'13207284'%IVP, %R_rho_S2T=2, tau=0.1, W0=1000 initial condition only 1 layer
                   };
               
        %This are new computation with the 512 grid point... 1024 is
        %difficult for eigenvalue computation, dedalus will be out of
        %memory.
        %Nz=512;
        slurm_num={'13207457',%R_rho_S2T=2, tau=0.1, W0=1 initial condition only 1 layer
                   '13207473'%R_rho_S2T=2, tau=0.1, W0=500.. (Note for different grid, the initial guess converge to staircase is different) initial condition, solution show the three layer structure
                   '13207459'%R_rho_S2T=2, tau=0.1, W0=1000, shows asymmetric solution...
                   };       
               
        %this is for Nz=128       
        slurm_num={'13208705',%one layer solution, R_rho_S2T=2, tau=0.1 + EVP
                   '13208739',%staircase 2 layer solution, R_rho_S2T=2, tau=0.1 + EVP
                   '13208744',%asymmetric solution, R_rho_S2T=2, tau=0.1 + EVP
                   '13207556',%IVP from 13207555
                   '13207561',%IVP from 13207553
                   '13207562',%IVP from 13207554
                   '13208655',%tau=0.01, R_rho_S2T=2, W0=1
                   '13208656',%tau=0.01,R_rho_S2T=2, W0=10
                   '13208657',%tau=0.01,R_rho_S2T=2, W0=100
                   '13208658',%tau=0.01,R_rho_S2T=2, W0=1000
                   '13208758',%tau=0.01,R_rho_S2T=2, 1-layer solution, continued from tau=0.1
                   '13208908',%IVP from 13208655
                   '13208911',%IVP from 13208656
                   '13208912',%IVP from 13208657
                   '13208914',%IVP from 13208658
                   '13208916',%IVP from 13208758
                   }; 
        %Nz=128, wavenumber from yang's scaling...  tau =0.01, Pr=7, Ra_T=100000, R_rho_T2S=2      
        slurm_num={'13209992',%continued from 13208655
                   '13209993',%continued from 13208656
                   '13209994',%continued from 13208657
                   '13209995',%continued from 13208658
                   '13209996',%continued from 13208758
                   '13210532',%another one layer kind of asymmetric solution, from initial guess
                   '13210528',%3 layer solution, continued from 13208658, but with kx as linspace(1, kx_final,10) so I can preserve symmetric...
                   '13209068',%kx continuation from 13208655
                   '13209069',%kx continuation from 13208656
                   '13209070',%kx continuation from 13208657
                   '13209071',%kx continuation from 13208658
                   '13209072' %kx continuation from 13208758
                     };  
%         slurm_num={'13209071'};
%         slurm_num=slurm_num(end-9:end-5);
%         slurm_num=slurm_num(end-13:end);

    case 'HB_benard_salt_finger_profile'
        slurm_num={'13209992',%continued from 13208655
                   '13209993',%continued from 13208656
                   '13209994',%continued from 13208657
                   '13209995',%continued from 13208658
                   '13209996',%continued from 13208758
                   '13210532',%another one layer kind of asymmetric solution, from initial guess
                   '13210528',%3 layer solution, continued from 13208658, but with kx as linspace(1, kx_final,10) so I can preserve symmetric...
                   };
    case 'HB_benard_salt_finger_profile_S1_yang'
        slurm_num={'13220808'};
    case 'HB_benard_salt_finger_kx'
        slurm_num={'13217936',%%kx increasing,[10:1:60], from 13209992, S2
                   %'13217937',%%kx increasing,[10:1:60], from 13209993, A2, branch jump, needs to be replcaed as 13220800
                   '13220800',%kx increasing, dense grid, from 13209993, A2,
                   '13217938',%%kx increasing,[10:1:60], from 13209994, S2L
                   %'13217939',%%kx increasing,[10:1:60], from 13209995, A3, branch jump, needs to be replaced as 13219318
                   '13219318',%kx increasing, 10 to 50, dense grid... from 13209995
                   '13217940',%%kx increasing,[10:1:60], from 13209996, S1
                   '13217941',%%kx increasing,[10:1:60], from 13210532, A1
                   '13217942',%%kx increasing,[10:1:60], from 13210528, S3
                   ...
                   '13217945',%%kx decreasing, [9.7815:0.1:20], from 13209992, S2
                   %'13217946',%%kx decreasing, [9.7815:0.1:20], from 13209993, A2, branch jump, replaced as 13219282
                   '13219282',%%kx decreasing, 9.7815 to 0.1, dense grid, A2 from 13209993
                   '13217947',%%kx decreasing, [9.7815:0.1:20], from 13209994, S2L
                   '13217948',%%kx decreasing, [9.7815:0.1:20], from 13209995, A3
                   '13217949',%%kx decreasing, [9.7815:0.1:20], from 13209996, S1
                   %'13217950',%%kx decreasing, [9.7815:0.1:20], from 13210532, A1, branch jump, replaced as 13219306
                   '13220804',%%kx decreasing, 9.7815 to 0.1, dense grid, A1 from 13210532
                   %'13217951',%%kx decreasing, [9.7815:0.1:20], from 13210528, S3, branch jump, replaced as 13219307
                   '13219307' %kx decreasing, 0.97815 to 0.1, dense grid
                   };
    case 'HB_benard_salt_finger_R_rho_T2S'    
        slurm_num={ %%next one is continuing with increasing R_rho_T2S
                   '13210696',%%R_rho_T2S increasing, from 13209992
                   '13210702',%%R_rho_T2S increasing, from 13209993
                   '13210703',%%R_rho_T2S increasing, from 13209994
                   '13210704',%%R_rho_T2S increasing, from 13209995
                   '13210705',%%R_rho_T2S increasing, from 13209996
                   '13210706',%%R_rho_T2S increasing, from 13210532
                   '13210707',%%R_rho_T2S increasing, from 13210528
                   ... %%next one is continuing with decreasing R_rho_T2S
                   '13210709',%%R_rho_T2S decreasing, from 13209992
                   '13210710',%%R_rho_T2S decreasing, from 13209993
                   '13210711',%%R_rho_T2S decreasing, from 13209994
                   '13210712',%%R_rho_T2S decreasing, from 13209995
                   '13210713',%%R_rho_T2S decreasing, from 13209996
                   '13210714',%%R_rho_T2S decreasing, from 13210532
                   '13210715',%%R_rho_T2S decreasing, from 13210528
                   ...
                   '13212083',%%R_rho_T2S increasing, from 10-500 from 13210696
                   '13212084',%%R_rho_T2S increasing, from 10-500 from 13210702
                   '13212085',%%R_rho_T2S increasing, from 10-500 from 13210703
                   '13212086',%%R_rho_T2S increasing, from 10-500 from 13210704
                   '13212087',%%R_rho_T2S increasing, from 10-500 from 13210705
                   '13212088',%%R_rho_T2S increasing, from 10-500 from 13210706
                   '13212089'%%R_rho_T2S increasing, from 10-500 from 13210707
                   }
%           slurm_num=slurm_num(end);
    case 'HB_benard_salt_finger_tau'
               slurm_num={'13212525',%tau continuation, from 13209992
                          '13212526',%tau continuation, from 13209993
                          '13212527',%tau continuation, from 13209994
                          '13212528',%tau continuation, from 13209995
                          '13212529',%tau continuation, from 13209996
                          '13212530',%tau continuation, from 13210532
                          '13212531' %tau continuation, from 13210528
                          };
                      
               %I need to check what they actually continue from
               slurm_num={'13217806',%tau continuation 0.01:0.01:1, from 13209992
                          '13217807',%tau continuation, from 13209993
                          '13217808',%tau continuation, from 13209994
                          '13217809',%tau continuation, from 13209995
                          '13217810',%tau continuation, from 13209996
                          '13217812',%tau continuation, from 13210532
                          '13217813',%tau continuation, from 13210528
                          };
               slurm_num=slurm_num(end);
                      
    case 'HB_benard_diffusive_kx'
        
        slurm_num={'13175327',
                    '13175339'};
        slurm_num=slurm_num(2);
    case 'HB_benard_diffusive_shear'
        slurm_num={'13114223'
                   };
    case 'hewitt_2_layer_Omega'
        slurm_num={'12829183'};
        %slurm_num={'12829121'}
    case 'herring_1963_free'
        slurm_num={'12761527',...Ra=4*10^3, kx=0.8pi
            '12761528',... Ra=10^4, kx=pi
            '12761529',...Ra=10^5, kx=1.5pi
            '12761536',...Ra=10^6, kx=1.5pi
            '12761538',...Ra=10^6, kx=6pi
            '12761539',...Ra=10^6, kx=9pi 
                    }
    case 'herring_1964_rigid'  
        %slurm_num={'12760848'}; continuation from this case, Ra_T=10^6,
        %kx=1
        slurm_num={'12761520',...Ra=4*10^3, kx=3
                '12761519',...Ra=10^4, kx=3
                '12761518',... Ra=10^5, kx=5
                '12761515'%Ra=10^6, kx=5.25
                }
    case 'trevisan_tau'
        slurm_num={'13037570',
                    '13037571',
                    '13037572',
                    '13037573',
                    '13037574'
                    };
    case 'trevisan_contour'
        slurm_num={'13039828'};
    case 'rosenberg_tau'
        slurm_num={'13060193',
                   '13060194',
                   '13060195',
                   '13060196'};
    case 'rosenberg_R_rho'
        slurm_num={'13060750',
                   '13060751',
                   '13060752',
                   '13060753'};
    case 'rosenberg_time_dependent'
        slurm_num={'13094886'};
    case 'mamou_contour_dirichlet'
        slurm_num={'13062713'};
    case 'mamou_contour_neumann'
        slurm_num={'13064536',
                   '13064537'};
    case 'mamou_contour_neumann'
    case 'test'
        slurm_num={'12834379'};
        %slurm_num={'12821151'};
        %slurm_num={'12760848'}; %Ra=10^6, kx=1
        
                %'12760763',
                %'12760764'};
        %slurm_num={'12760848',%for a range of Ra from 10^4 to 10^6
        %        '12760884'};%for Ra=10^6, but wavenumber from 1 to 20
        
    case 'HB_porous_mamou_Ra_T_55_SW'
        %Ra_T=55, Ra_S2T=5.5, Le=5, single-mode DNS of porous medium
        %oscillatory double-diffusive convection
%         slurm_num={'14446324'};
        slurm_num={'14493285'};
    case 'HB_porous_mamou_Ra_T_53_SW'
        %Ra_T=53, Le=5, single-mode DNS of oscillatory double-diffusive
        %convection
        %'14501184',...Ra_S2T=3.6924243
        %'14501185',...Ra_S2T=3.6924247
        
        %'14502039',...Ra_S2T=3.692424135
        slurm_num={'14502056',...Ra_S2T=3.6924241348
            '14502037',...Ra_S2T=3.69242414
            '14502032',...Ra_S2T=3.69242415
            '14501186',...Ra_S2T=3.6924242
            '14501183',...Ra_S2T=3.6924245
            '14501175',...Ra_S2T=3.692425
            '14501178',...Ra_S2T=3.692427
            '14501168',...Ra_S2T=3.69243
            '14501167',...Ra_S2T=3.69244
            '14500516',...Ra_S2T=3.69245
            '14500129',...Ra_S2T=3.6925
            '14500128',...Ra_S2T=3.6926
            '14500127',...Ra_S2T=3.6928
            '14500126',...Ra_S2T=3.693
            '14500125',...Ra_S2T=3.694
            '14493257',...Ra_S2T=3.695
            '14450520',...Ra_S2T=3.7
            '14450521',...Ra_S2T=3.8
            '14450507',...Ra_S2T=3.9
            '14451373',...Ra_S2T=4
            '14451370',...Ra_S2T=4.1
            '14450502',...Ra_S2T=4.2
            '14451369',...Ra_S2T=4.3
            '14451368',...Ra_S2T=4.4
            '14450501',...Ra_S2T=4.5
            '14450500',...Ra_S2T=4.6
            '14450499',...Ra_S2T=4.7
            '14450477',...Ra_S2T=4.8
            '14450478',...Ra_S2T=4.9
            '14450479',...Ra_S2T=5
            '14450481',...Ra_S2T=5.1
            '14450482',...Ra_S2T=5.2
            '14450483',...Ra_S2T=5.3
            '14450484',...Ra_S2T=5.4
            '14450506',...Ra_S2T=5.5
            '14450497'...Ra_S2T=5.6
            };
%             '14450487',...Ra_S2T=5.7
%             '14450489',...Ra_S2T=5.8
%             '14450491',...Ra_S2T=5.9
%             '14450492'...Ra_S2T=6.0
    case 'HB_porous_rosenberg_Ra_T_100_SW'
        % '14502705',...Ra_S2T=40.52
%         '14502728',...Ra_S2T=40.512
        % '14502729',...Ra_S2T=40.511
    %            '14502719',...Ra_S2T=40.515
%            '14502705',...Ra_S2T=40.52
        slurm_num={'14502744',...Ra_S2T=40.5105
            '14502726',...Ra_S2T=40.513
            '14502727',...Ra_S2T=40.517
            '14502659',...Ra_S2T=40.53
            '14502658',...Ra_S2T=40.55
            '14502656',...Ra_S2T=40.6
            '14502655',...Ra_S2T=40.7
            '14451123',...Ra_S2T=41
            '14451121',...Ra_S2T=43
            '14451120',...Ra_S2T=45
            '14451095',...Ra_S2T=50
            '14451094'...Ra_S2T=55
            };
        %'14451105',...Ra_S2T=5
        %'14451104',...Ra_S2T=10
        %'14451103',...Ra_S2T=15
        %'14451102',...Ra_S2T=20
        %'14451101',...Ra_S2T=25
        %'14451100',...Ra_S2T=30
end

%Update 2022/05/04, single mode DNS
% slurm_num={'13631229',...
%             '13631083',...
%             '13631099'};
% % slurm_num={'13631270',...
% %             '13631271',...
% %             '13631272'};
% slurm_num={'13632918',...
%             '13632920',...
%             '13632921'};
% 
% %single mode DNS with initial condition as the traveling wave        
% slurm_num={'13633162',...: kx=3
%             '13633163',...: kx=4
%             '13633164',...: kx=5
%             '13633165',...: kx=6
%             '13633176',...: kx=7 
%             '13633177',...: kx=8
%             '13633178',...: kx=9
%             '13633179',...: kx=10
%             '13633180',...: kx=11
%             '13633181',...: kx=12
%             '13633182',...: kx=13
%             '13633191',...: kx=14
%             '13633194',...: kx=15
%             '13633198'...: kx=18
%             };
% 
% % %single mode DNS with initial condition as the tilted finger
% % slurm_num={'13633210',...: kx=2
% %         '13633211',...: kx=3
% %         '13633212',...: kx=4
% %         '13633213',...: kx=5
% %         '13633214',...: kx=6
% %         '13633215',...: kx=7
% %         '13633216',...: kx=8
% %         '13633218',...: kx=9
% %         '13633220',...: kx=10
% %         '13633221',...: kx=11
% %         '13633222',...: kx=12
% %         '13633223',...: kx=13
% %         '13633224',...: kx=14
% %         '13633225',...: kx=15
% %         '13633226',...: kx=16
% %         '13633227'...: kx=17
% %         };
% %         
% slurm_num={'13633381',...
%             '13633382'};
% slurm_num={'13633224'};


%standing wave at Ra_T=53..

% group_name=[];
flag.print=0;
flag.visible=0;
flag.video=0;
flag.no_ylabel=0;
flag.post_plot=1;
for slurm_ind=1:length(slurm_num)
    content=dir(['C:\Data\dedalus\dedalus_',...
        slurm_num{slurm_ind}]);
    content(1:3)=[];%remove these...
    for content_ind=1:length(content)
        content_ind
        if strcmp(content(content_ind).name(1:8),'analysis')
            h5_name=['C:\Data\dedalus\dedalus_',...
            slurm_num{slurm_ind},...
                '\',content(content_ind).name,'\analysis_s1.h5'];

            dedalus_post_my{slurm_ind,content_ind}=dedalus_post(h5_name,flag);
            dedalus_post_my{slurm_ind,content_ind}.uvw_hewitt=0;
            %data_Nu{1}.x(slurm_ind,content_ind)=dedalus_post_my{slurm_ind,content_ind}.kx;
            %data_Nu{1}.y(slurm_ind,content_ind)=dedalus_post_my{slurm_ind,content_ind}.Ra_T;
            %data_Nu{1}.z(slurm_ind,content_ind)=dedalus_post_my{slurm_ind,content_ind}.Nu(1);
            switch dedalus_post_my{slurm_ind,content_ind}.problem
                case {'BVP','IVP','EVP'}
                    dedalus_post_my{slurm_ind,content_ind}=dedalus_post_my{slurm_ind,content_ind}.dedalus_post_bvp();
                    if flag.print || flag.visible
                        dedalus_post_my{slurm_ind,content_ind}=dedalus_post_my{slurm_ind,content_ind}.bvp_plot;
                    end
%                 case 'EVP'
%                     dedalus_post_my{slurm_ind,content_ind}=dedalus_post_my{slurm_ind,content_ind}.dedalus_post_evp();
%                     if flag.print || flag.visible
%                         dedalus_post_my{slurm_ind,content_ind}=dedalus_post_my{slurm_ind,content_ind}.evp_plot;
%                     end
%                 case 'IVP'
                otherwise
                    error('Wrong problem for dedalus');
            end
        end
    end
end

%%------old plotting
% if flag.post_plot
% 
%     if dedalus_post_my{slurm_ind,content_ind}.dy_T_mean<0
%         background_T=1-dedalus_post_my{slurm_ind,content_ind}.z_list;
%     elseif dedalus_post_my{slurm_ind,content_ind}.dy_T_mean>0
%         background_T=dedalus_post_my{slurm_ind,content_ind}.z_list;
%     end
%     switch group_name
% %                     case 'hewitt_2D'
% %                         data_z{1}.x(slurm_ind)=dedalus_post_my{slurm_ind}.Ra_T;
% %                         data_z{1}.y(slurm_ind)=dedalus_post_my{slurm_ind}.z_T_BL;
% %                         data_z{2}.x(slurm_ind)=dedalus_post_my{slurm_ind}.Ra_T;
% %                         data_z{2}.y(slurm_ind)=dedalus_post_my{slurm_ind}.z_T_rms_max;
% %                         %data_T_BL{slurm_ind}.x=deda
% %                         %data_S_BL{slurm_ind}.x=dedalus_post_my{slurm_ind}.Ra_T;
% %                         %data_S_BL{slurm_ind}.y=dedalus_post_my{slurm_ind}.S_BL;
% %                     case 'hewitt_3D'
% %                         
%         case {'HB_porous_diffusive_BC','HB_porous_finger_BC'}
%           if dedalus_post_my{slurm_ind,content_ind}.dy_S_mean<0
%               background_S=1-dedalus_post_my{slurm_ind,content_ind}.z_list;
%           elseif dedalus_post_my{slurm_ind,content_ind}.dy_S_mean>0
%               background_S=dedalus_post_my{slurm_ind,content_ind}.z_list;
%           end
%           data_S{slurm_ind,content_ind}.y=dedalus_post_my{slurm_ind,content_ind}.z_list;
%           data_S{slurm_ind,content_ind}.x=dedalus_post_my{slurm_ind,content_ind}.S_0+background_S;
%           data_S{slurm_ind,content_ind}.x=(data_S{slurm_ind,content_ind}.x-min(data_S{slurm_ind,content_ind}.x))/(max(data_S{slurm_ind,content_ind}.x)-min(data_S{slurm_ind,content_ind}.x));            
%         case 'HB_benard_kx'
%             data_Nu{1}.x(content_ind)=dedalus_post_my{slurm_ind,content_ind}.kx;
%             data_Nu{1}.y(content_ind)=dedalus_post_my{slurm_ind,content_ind}.Nu(1);
%         case 'HB_benard_Ra'
%             data_Nu{1}.x(content_ind)=dedalus_post_my{slurm_ind,content_ind}.Ra_T;
%             data_Nu{1}.y(content_ind)=dedalus_post_my{slurm_ind,content_ind}.Nu(1);
%         case 'HB_porous_kx'
%             data_Nu{1}.x(content_ind)=dedalus_post_my{slurm_ind,content_ind}.kx;
%             data_Nu{1}.y(content_ind)=dedalus_post_my{slurm_ind,content_ind}.Nu(1);
%         case 'HB_porous_Nu_kx_Ra'
%             %data_Nu{1}.x()=dedalus_post_my{slurm_ind,content_ind}.kx;
%             %data_Nu{1}.y(slurm_ind)=dedalus_post_my{slurm_ind,content_ind}.Ra_T;
%             data_Nu{1}.z(slurm_ind,content_ind)=dedalus_post_my{slurm_ind,content_ind}.Nu(1);
%     end
% end

%%----


%%plotting
if flag.post_plot
    switch group_name
        case {'HB_porous_mamou_Ra_T_53_SW',...
                'HB_porous_rosenberg_Ra_T_100_SW'}
            for slurm_ind=1:length(slurm_num)
                Ra_S2T_list(slurm_ind,1)=dedalus_post_my{slurm_ind}.Ra_S2T;
                Nz=dedalus_post_my{slurm_ind}.Nz;
                w_hat_diff=diff(dedalus_post_my{slurm_ind}.w_hat(Nz/2,:));
                ind_change=find(diff(sign(w_hat_diff)));
                ind_later=2*round(length(ind_change)/4);
                Nu_mean(slurm_ind,1)=mean(dedalus_post_my{slurm_ind}.Nu(1,ind_change(ind_later):ind_change(end)));
                Sh_mean(slurm_ind,1)=mean(dedalus_post_my{slurm_ind}.Nu_S(1,ind_change(ind_later):ind_change(end)));
                period(slurm_ind,1)=(dedalus_post_my{slurm_ind}.t_list(ind_change(end))-dedalus_post_my{slurm_ind}.t_list(ind_change(ind_later)))/((length(ind_change)-ind_later)/2);
                
                %new method using fft to get the oscillation frequency
                t_ind_begin=length(dedalus_post_my{slurm_ind}.t_list)/4;
                t_list_dedalus=dedalus_post_my{slurm_ind}.t_list(t_ind_begin:end);
                dt=mean(diff(dedalus_post_my{slurm_ind}.t_list(t_ind_begin:end)));
                Nt=length(t_list_dedalus);
                Fs=1/dt;
                freq=Fs*(0:(Nt/2))/Nt;
                spec_tmp=abs(fft(dedalus_post_my{slurm_ind}.w_hat(Nz/2,t_ind_begin:end))/Nt);
                spec_tmp=spec_tmp(1:Nt/2+1);
                spec_tmp(2:end-1)=2*spec_tmp(2:end-1);
                [val,ind]=max(spec_tmp);
                Nu_mean(slurm_ind,1)=mean(dedalus_post_my{slurm_ind}.Nu(1,t_ind_begin:end));
                Sh_mean(slurm_ind,1)=mean(dedalus_post_my{slurm_ind}.Nu_S(1,t_ind_begin:end));
                period(slurm_ind,1)=1/freq(ind);
            end
            
        case 'HB_porous_mamou_Ra_T_55_SW'
            %post-processing for the standing wave.
            %Update 2022/10/07
            folder_name=['C:\Data\dedalus\dedalus_',...
                        slurm_num{slurm_ind},'/analysis/'];

%             ind_start=1306;
%             ind_end=1463; %slurm_num=14446324
%             
            ind_start=13076;
            ind_end=14644; %slurm_num=14493285
            
            kx=dedalus_post_my{1}.kx;
            phi_hat=dedalus_post_my{1}.w_hat/(1i*kx);
            data{1}.x=dedalus_post_my{1}.t_list(ind_start:ind_end);
            data{1}.y=10*(dedalus_post_my{1}.Nu(1,ind_start:ind_end)-1);
            data{2}.x=dedalus_post_my{1}.t_list(ind_start:ind_end);
            data{2}.y=dedalus_post_my{1}.Nu_S(1,ind_start:ind_end);
            data{3}.x=dedalus_post_my{1}.t_list(ind_start:ind_end);
            data{3}.y=2*(imag(phi_hat(64,ind_start:ind_end)));
            plot_config.legend_list={1,'$10[nu(t)-1]$','$sh(t)$', '$\psi_{\mathrm {mid}}(t)$'};
            plot_config.name=[folder_name,group_name,'nu_sh_phi_m.png'];
            plot_config.label_list={1,'$t$',''};
            plot_config.print=1;
            plot_config.fontsize_legend=28;
            plot_config.linewidth=3;
            plot_line(data,plot_config);
            
            ind_interval=ind_end-ind_start;
            streamfunction_ind_list=[ind_start,...
                    ind_start+round(0.2697*ind_interval),...
                    ind_start+round(0.2714*ind_interval),...
                    ind_start+round(0.27166*ind_interval),...
                    ind_start+round(0.2723*ind_interval),...
                    ind_start+round(0.5*ind_interval)];
       
            x=linspace(0,2*pi,1000);
            y=dedalus_post_my{1}.z_list;
            for streamfunction_ind=streamfunction_ind_list

                data{2}.x=x; data{2}.y=y;
                data{2}.z=2*real(dedalus_post_my{1}.w_hat(:,streamfunction_ind)/(1i*dedalus_post_my{1}.kx)*exp(1i*x));
                [data{1}.x,data{1}.y]=meshgrid(x,y);
                data{1}.z=NaN*ones(size(data{1}.x));

                plot_config.xlim_list=[1,0,2*pi];
                plot_config.xtick_list=[1,0,pi/2,pi,3*pi/2,2*pi];
                plot_config.xticklabels_list={1,'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'};
                plot_config.ylim_list=[1,0,1];
                plot_config.label_list={1,'$x k_x$','$z$'};
                plot_config.streamline=2;
                plot_config.user_color_style_marker_list={'k-','r--'};
                plot_config.panel_num=2;
                plot_config.arrow_ratio=0.8;
                plot_config.linewidth=3;
                plot_config.colorbar=0;
                plot_config.fontsize=28;
                plot_config.print_size=[1,500,900];
                plot_config.name=[folder_name,group_name,'streamfunction_t_ind=',num2str(streamfunction_ind),'.png'];
                plot_contour(data,plot_config);
            end
            
        case 'hewitt_2D'
            %This is plotting the validation of comparison with DNS mean
            %temperature at three different Rayleigh
            %number available from Hewitt et al. (2012)
            for slurm_ind=1:size(dedalus_post_my,1)
                for content_ind=1:size(dedalus_post_my,2)
                    data_T{slurm_ind,content_ind}.y=dedalus_post_my{slurm_ind,content_ind}.z_list;
                    data_T{slurm_ind,content_ind}.x=dedalus_post_my{slurm_ind,content_ind}.T_0+background_T;
                    data_T{slurm_ind,content_ind}.x=(data_T{slurm_ind,content_ind}.x-min(data_T{slurm_ind,content_ind}.x))/(max(data_T{slurm_ind,content_ind}.x)-min(data_T{slurm_ind,content_ind}.x));
                end
            end
            porous_hewitt_2D=dedalus_post_my{1}.get_porous_hewitt_2D();
            sparse_ind=5;
            data_T{4}.x=porous_hewitt_2D.Ra_10000_z_T_0(1:sparse_ind:length(porous_hewitt_2D.Ra_10000_z_T_0),2);
            data_T{4}.y=porous_hewitt_2D.Ra_10000_z_T_0(1:sparse_ind:length(porous_hewitt_2D.Ra_10000_z_T_0),1);
            data_T{5}.x=porous_hewitt_2D.Ra_20000_z_T_0(1:sparse_ind:length(porous_hewitt_2D.Ra_20000_z_T_0),2);
            data_T{5}.y=porous_hewitt_2D.Ra_20000_z_T_0(1:sparse_ind:length(porous_hewitt_2D.Ra_20000_z_T_0),1);
            data_T{6}.x=porous_hewitt_2D.Ra_40000_z_T_0(1:sparse_ind:length(porous_hewitt_2D.Ra_40000_z_T_0),2);
            data_T{6}.y=porous_hewitt_2D.Ra_40000_z_T_0(1:sparse_ind:length(porous_hewitt_2D.Ra_40000_z_T_0),1);

            sparse_ind=1;
            data_T{7}.x=porous_hewitt_2D.Ra_10000_z_T_0(1:sparse_ind:length(porous_hewitt_2D.Ra_10000_z_T_0),2);
            data_T{7}.y=porous_hewitt_2D.Ra_10000_z_T_0(1:sparse_ind:length(porous_hewitt_2D.Ra_10000_z_T_0),1);
            data_T{8}.x=porous_hewitt_2D.Ra_20000_z_T_0(1:sparse_ind:length(porous_hewitt_2D.Ra_20000_z_T_0),2);
            data_T{8}.y=porous_hewitt_2D.Ra_20000_z_T_0(1:sparse_ind:length(porous_hewitt_2D.Ra_20000_z_T_0),1);
            data_T{9}.x=porous_hewitt_2D.Ra_40000_z_T_0(1:sparse_ind:length(porous_hewitt_2D.Ra_40000_z_T_0),2);
            data_T{9}.y=porous_hewitt_2D.Ra_40000_z_T_0(1:sparse_ind:length(porous_hewitt_2D.Ra_40000_z_T_0),1);

            %This is not used but just for the legend plotting
            data_T{10}.x=NaN;
            data_T{10}.y=NaN;
            data_T{11}.x=NaN;
            data_T{11}.y=NaN;
            data_T{12}.x=NaN;
            data_T{12}.y=NaN;
            
            %plot_config.legend_list={1,'$Ra_T=10000$','$Ra_T=20000$','$Ra_T=40000$'};
            plot_config.label_list={1,'$\bar{T}_0+1-z$','$z$'};
            plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'T_0.png'];
            plot_config.legend_list={1,'Ra=10000','Ra=20000', 'Ra=40000',...
                    'Ra=10000 (DNS)', 'Ra=20000 (DNS)','Ra=40000 (DNS)'};
            plot_config.legend_index=[1,2,3,10,11,12];
            plot_config.fontsize_legend=26;
            plot_config.fontsize=32;
            plot_config.linwidth=3;
            plot_config.Markerindex=3;
            plot_config.user_color_style_marker_list={'r-','k--','b-.','rsquare','k^','bo','r-','k--','b-.','r-square','k--^','b-.o'};
            plot_config.print_size=[1,900,1000];
            plot_config.linwidth=3;
%             plot_config.legend_list={0};
            plot_config.xlim_list=[1,0.45,0.55];
            plot_line(data_T,plot_config);
%plot_config.legend_list={1,'$Ra=10^4$','$Ra=2\times 10^4$','$Ra_T=4\times 10^4$','$Ra=10^4$ Hewitt et al. (2012)','$Ra=2\times 10^4$ Hewitt et al. (2012)','$Ra_T=4\times 10^4$ Hewitt et al. (2012)' };

%             plot_config.user_color_style_marker_list={'r-','k--','b-.','rsquare','k^','bo'};
%             plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'T_0_local.png'];
%             plot_line(data_T,plot_config);
%             
% 
%             plot_config.loglog=[1,1];
%             plot_config.xlim_list=0; 
%             plot_config.ylim_list=0
%             plot_config.label_list={1,'$Ra_T$','$z_T$'};
%             [eta,c0]=scaling(data_z{1}.x,data_z{1}.y);
%             plot_line(data_z,plot_config);

        case 'hewitt_3D'
            %This is plotting the validation of comparison with DNS mean
            %temperature and the rms value at three different Rayleigh
            %number available from Hewitt et al. (2012)
            for slurm_ind=1:size(dedalus_post_my,1)
                for content_ind=1:size(dedalus_post_my,2)
                    data_T{slurm_ind,content_ind}.y=dedalus_post_my{slurm_ind,content_ind}.z_list;
                    data_T{slurm_ind,content_ind}.x=dedalus_post_my{slurm_ind,content_ind}.T_0+background_T;
                    data_T{slurm_ind,content_ind}.x=(data_T{slurm_ind,content_ind}.x-min(data_T{slurm_ind,content_ind}.x))/(max(data_T{slurm_ind,content_ind}.x)-min(data_T{slurm_ind,content_ind}.x));
                end
            end
            porous_hewitt_3D=dedalus_post_my{1}.get_porous_hewitt_3D();
            sparse_ind=5;
            data_T{4}.x=porous_hewitt_3D.Ra_4000_z_T_0(1:sparse_ind:length(porous_hewitt_3D.Ra_4000_z_T_0),2);
            data_T{4}.y=porous_hewitt_3D.Ra_4000_z_T_0(1:sparse_ind:length(porous_hewitt_3D.Ra_4000_z_T_0),1);
            data_T{5}.x=porous_hewitt_3D.Ra_8000_z_T_0(1:sparse_ind:length(porous_hewitt_3D.Ra_8000_z_T_0),2);
            data_T{5}.y=porous_hewitt_3D.Ra_8000_z_T_0(1:sparse_ind:length(porous_hewitt_3D.Ra_8000_z_T_0),1);
            data_T{6}.x=porous_hewitt_3D.Ra_16000_z_T_0(1:sparse_ind:length(porous_hewitt_3D.Ra_16000_z_T_0),2);
            data_T{6}.y=porous_hewitt_3D.Ra_16000_z_T_0(1:sparse_ind:length(porous_hewitt_3D.Ra_16000_z_T_0),1);

            sparse_ind=1;
            data_T{7}.x=porous_hewitt_3D.Ra_4000_z_T_0(1:sparse_ind:length(porous_hewitt_3D.Ra_4000_z_T_0),2);
            data_T{7}.y=porous_hewitt_3D.Ra_4000_z_T_0(1:sparse_ind:length(porous_hewitt_3D.Ra_4000_z_T_0),1);
            data_T{8}.x=porous_hewitt_3D.Ra_8000_z_T_0(1:sparse_ind:length(porous_hewitt_3D.Ra_8000_z_T_0),2);
            data_T{8}.y=porous_hewitt_3D.Ra_8000_z_T_0(1:sparse_ind:length(porous_hewitt_3D.Ra_8000_z_T_0),1);
            data_T{9}.x=porous_hewitt_3D.Ra_16000_z_T_0(1:sparse_ind:length(porous_hewitt_3D.Ra_16000_z_T_0),2);
            data_T{9}.y=porous_hewitt_3D.Ra_16000_z_T_0(1:sparse_ind:length(porous_hewitt_3D.Ra_16000_z_T_0),1);

            %
            data_T{10}.x=NaN;
            data_T{10}.y=NaN;
            data_T{11}.x=NaN;
            data_T{11}.y=NaN;
            data_T{12}.x=NaN;
            data_T{12}.y=NaN;
            %plot_config.legend_list={1,'$Ra_T=4000$','$Ra_T=8000$','$Ra_T=16000$'};
            plot_config.label_list={1,'$\bar{T}_0+1-z$','$z$'};
            plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'T_0.png'];
            plot_config.legend_list={1,'Ra=4000','Ra=8000', 'Ra=16000',...
                    'Ra=4000 (DNS)', 'Ra=8000 (DNS)','Ra=16000 (DNS)'};
            plot_config.legend_index=[1,2,3,10,11,12];
            plot_config.fontsize_legend=26;
            plot_config.fontsize=32;
            plot_config.linwidth=3;
            plot_config.Markerindex=3;
            plot_config.user_color_style_marker_list={'r-','k--','b-.','rsquare','k^','bo','r-','k--','b-.','r-square','k--^','b-.o'};
            plot_config.print_size=[1,900,1000];
            plot_config.xlim_list=[1,0.45,0.55];
            plot_line(data_T,plot_config);
            
            %plot of T_rms comparison...
            rms_list={'T','w','u'};
            clear data_T
            plot_config.linwidth=3;
            plot_config.fontsize=32;
            plot_config.legend_list={0};
            for rms_ind=1:length(rms_list)
                rms_name=rms_list{rms_ind};
                if rms_name=='u'
                    suffix='_tilde';
                else 
                    suffix='_hat';
                end
                data_T{1}.x=abs(dedalus_post_my{1}.([rms_name,suffix]))*sqrt(2);
                data_T{1}.y=dedalus_post_my{1}.z_list;
                data_T{2}.x=abs(dedalus_post_my{2}.([rms_name,suffix]))*sqrt(2);
                data_T{2}.y=dedalus_post_my{2}.z_list;
                data_T{3}.x=abs(dedalus_post_my{3}.([rms_name,suffix]))*sqrt(2);
                data_T{3}.y=dedalus_post_my{3}.z_list;

                porous_hewitt_3D=dedalus_post_my{1}.get_porous_hewitt_3D();
                sparse_ind=1;
                data_T{4}.x=porous_hewitt_3D.(['Ra_4000_z_',rms_name,'_rms'])(1:sparse_ind:length(porous_hewitt_3D.(['Ra_4000_z_',rms_name,'_rms'])),2);
                data_T{4}.y=porous_hewitt_3D.(['Ra_4000_z_',rms_name,'_rms'])(1:sparse_ind:length(porous_hewitt_3D.(['Ra_4000_z_',rms_name,'_rms'])),1);
                data_T{5}.x=porous_hewitt_3D.(['Ra_8000_z_',rms_name,'_rms'])(1:sparse_ind:length(porous_hewitt_3D.(['Ra_8000_z_',rms_name,'_rms'])),2);
                data_T{5}.y=porous_hewitt_3D.(['Ra_8000_z_',rms_name,'_rms'])(1:sparse_ind:length(porous_hewitt_3D.(['Ra_8000_z_',rms_name,'_rms'])),1);
                data_T{6}.x=porous_hewitt_3D.(['Ra_16000_z_',rms_name,'_rms'])(1:sparse_ind:length(porous_hewitt_3D.(['Ra_16000_z_',rms_name,'_rms'])),2);
                data_T{6}.y=porous_hewitt_3D.(['Ra_16000_z_',rms_name,'_rms'])(1:sparse_ind:length(porous_hewitt_3D.(['Ra_16000_z_',rms_name,'_rms'])),1);

                sparse_ind=5;
                data_T{7}.x=porous_hewitt_3D.(['Ra_4000_z_',rms_name,'_rms'])(1:sparse_ind:length(porous_hewitt_3D.(['Ra_4000_z_',rms_name,'_rms'])),2);
                data_T{7}.y=porous_hewitt_3D.(['Ra_4000_z_',rms_name,'_rms'])(1:sparse_ind:length(porous_hewitt_3D.(['Ra_4000_z_',rms_name,'_rms'])),1);
                data_T{8}.x=porous_hewitt_3D.(['Ra_8000_z_',rms_name,'_rms'])(1:sparse_ind:length(porous_hewitt_3D.(['Ra_8000_z_',rms_name,'_rms'])),2);
                data_T{8}.y=porous_hewitt_3D.(['Ra_8000_z_',rms_name,'_rms'])(1:sparse_ind:length(porous_hewitt_3D.(['Ra_8000_z_',rms_name,'_rms'])),1);
                data_T{9}.x=porous_hewitt_3D.(['Ra_16000_z_',rms_name,'_rms'])(1:sparse_ind:length(porous_hewitt_3D.(['Ra_16000_z_',rms_name,'_rms'])),2);
                data_T{9}.y=porous_hewitt_3D.(['Ra_16000_z_',rms_name,'_rms'])(1:sparse_ind:length(porous_hewitt_3D.(['Ra_16000_z_',rms_name,'_rms'])),1);

                plot_config.label_list={1,['$',rms_name,'_{\textrm{rms}}(z)$'],'$z$'};
                plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,rms_name,'_rms.png'];
                plot_config.Markerindex=3;
                plot_config.user_color_style_marker_list={'r-','k--','b-.','r-','k--','b-.','rsquare','k^','bo'};
                plot_config.print_size=[1,600,1000];
                plot_config.xlim_list=[0,0.45,0.55];
                plot_line(data_T,plot_config);
            end
            
            %plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'T_0_local.png'];
            %plot_line(data_T,plot_config);
        case 'hewitt_2D_Ra'
            %This is plotting the trend over Rayleigh number for 2D
            
            clear data_Nu
            %read the Rayleigh number and order the data by Rayleigh number
            for content_ind=1:length(dedalus_post_my)
                data_Nu{1}.x(content_ind)=dedalus_post_my{1,content_ind}.Ra_T;
            end
            [data_Nu{1}.x,ind]=sort(data_Nu{1}.x);
            dedalus_post_my=dedalus_post_my(1,ind);
            for content_ind=1:length(dedalus_post_my)
                data_Nu{1}.y(content_ind)=dedalus_post_my{1,content_ind}.Nu(1);
            end
            porous_hewitt_2D=dedalus_post_my{1}.get_porous_hewitt_2D();
            porous_wen_2D=dedalus_post_my{1}.get_porous_wen_2D();
            porous_otero_bound=dedalus_post_my{1}.get_porous_otero_bound();
            
            sparse_ind=4;
            data_Nu{2}.x=porous_hewitt_2D.Ra_Nu(1:sparse_ind:length(porous_hewitt_2D.Ra_Nu),1);
            data_Nu{2}.y=porous_hewitt_2D.Ra_Nu(1:sparse_ind:length(porous_hewitt_2D.Ra_Nu),2);
            data_Nu{3}.x=porous_wen_2D.Ra_Nu_DNS(1:sparse_ind:length(porous_wen_2D.Ra_Nu_DNS),1);
            data_Nu{3}.y=porous_wen_2D.Ra_Nu_DNS(1:sparse_ind:length(porous_wen_2D.Ra_Nu_DNS),2);
            
            sparse_ind=1;
            data_Nu{6}.x=porous_hewitt_2D.Ra_Nu(:,1);
            data_Nu{6}.y=porous_hewitt_2D.Ra_Nu(:,2);
            data_Nu{7}.x=porous_wen_2D.Ra_Nu_DNS(:,1);
            data_Nu{7}.y=porous_wen_2D.Ra_Nu_DNS(:,2);
            
            data_Nu{4}.x=porous_wen_2D.Ra_Nu_steady(:,1);
            data_Nu{4}.y=porous_wen_2D.Ra_Nu_steady(:,2);
            data_Nu{5}.x=porous_otero_bound.Ra_Nu_bound(:,1);
            data_Nu{5}.y=porous_otero_bound.Ra_Nu_bound(:,2);
            
            %For legend plotting
%             data_Nu{8}.x=NaN;
%             data_Nu{8}.y=NaN;
%             data_Nu{9}.x=NaN;
%             data_Nu{9}.y=NaN;
%             
            plot_config.linwidth=3;
            plot_config.label_list={1,'$Ra$','$Nu$'};
            plot_config.loglog=[1,1];
            plot_config.Markerindex=3;
            plot_config.xlim_list=0; plot_config.ylim_list=0;
            plot_config.print_size=[1,1000,1000];
            plot_config.user_color_style_marker_list={'k-','r^','rsquare','bo','b--','r-^','r-^'}; 
            plot_config.legend_list={1,'Single mode','DNS (Hewitt \it{et al.} \rm 2012)', 'DNS (Wen \it{et al.} \rm 2015)',...
               'Steady solutions (Wen \it{et al.} \rm 2015)','Upper bound (Otero \it{et al.} \rm 2004)'};
            plot_config.fontsize_legend=20;
            plot_config.fontsize=32;
            plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'Ra_Nu.png'];
            plot_config.linewidth=3;
            plot_config.xtick_list=[1,1,10,100,1000,10000,100000,1000000];
            plot_config.ytick_list=[1,1,10,100,1000];
            plot_config.ylim_list=[1,1,2500];
            plot_config.xlim_list=[1,10,1000000];
            plot_line(data_Nu,plot_config);
            [data_Nu_scaling.eta,data_Nu_scaling.c0]=scaling(data_Nu{1}.x(50:end),data_Nu{1}.y(50:end));

            %--------read the rms value and plot
            for content_ind=1:length(dedalus_post_my)
                data_rms{1}.y(content_ind)=dedalus_post_my{1,content_ind}.T_rms_mid;
                data_rms{2}.y(content_ind)=dedalus_post_my{1,content_ind}.w_rms_mid;
                data_rms{3}.y(content_ind)=dedalus_post_my{1,content_ind}.u_rms_mid;
            end
            
            data_rms{1}.x=data_Nu{1}.x;
            data_rms{2}.x=data_Nu{1}.x;
            data_rms{3}.x=data_Nu{1}.x;
            data_rms{4}.x=porous_hewitt_2D.Ra_T_rms(:,1);
            data_rms{4}.y=porous_hewitt_2D.Ra_T_rms(:,2);
            data_rms{5}.x=porous_hewitt_2D.Ra_w_rms(:,1);
            data_rms{5}.y=porous_hewitt_2D.Ra_w_rms(:,2);
            data_rms{6}.x=porous_hewitt_2D.Ra_u_rms(:,1);
            data_rms{6}.y=porous_hewitt_2D.Ra_u_rms(:,2);
            
            plot_config.label_list={1,'$Ra$',[]};
            plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'Ra_rms.png'];
            plot_config.Markerindex=3;
            plot_config.user_color_style_marker_list={'r-','b--','k-.','rsquare','bo','k^'};    
            plot_config.loglog=[1,0];
            plot_config.linewidth=3;
            plot_config.legend_list={1,'$T_{\textrm{rms}}(z=0.5)$','$w_{\textrm{rms}}(z=0.5)$','$u_{\textrm{rms}}(z=0.5)$',...
                '$T_{\textrm{rms}}(z=0.5)$ (DNS)','$w_{\textrm{rms}}(z=0.5)$ (DNS)','$u_{\textrm{rms}}(z=0.5)$ (DNS)'};
            plot_config.fontsize_legend=18;
            plot_config.fontsize=32;
            plot_config.ylim_list=[1,0,0.15];
            plot_config.xlim_list=[1,100,1000000];
            plot_config.ytick_list=[1,0,0.05,0.1,0.15];
            plot_line(data_rms,plot_config);
            
            for content_ind=1:length(dedalus_post_my)
                d_T_total=dedalus_post_my{1,content_ind}.d_T_0+dedalus_post_my{1,content_ind}.dy_T_mean;
                data_d_T{1}.y(content_ind)=dedalus_post_my{1,content_ind}.Nu(1);
                data_d_T{2}.y(content_ind)=max(d_T_total);
                data_d_T{3}.y(content_ind)=dedalus_post_my{1,content_ind}.T_rms_max;
                data_d_T{4}.y(content_ind)=-dedalus_post_my{1,content_ind}.d_T_0_mid;
                data_d_T{5}.y(content_ind)=dedalus_post_my{1,content_ind}.z_T_BL;
                data_d_T{6}.y(content_ind)=dedalus_post_my{1,content_ind}.z_T_rms_max;
                if data_d_T{6}.y(content_ind)>0.5
                    data_d_T{6}.y(content_ind)=1-data_d_T{6}.y(content_ind);
                end
            end
            negative_max_d_T_0=find(data_d_T{2}.y<0);
            data_d_T{5}.y(negative_max_d_T_0)=NaN;
            data_d_T{1}.x=data_Nu{1}.x;
            data_d_T{2}.x=data_Nu{1}.x;
            data_d_T{3}.x=data_Nu{1}.x;
            data_d_T{4}.x=data_Nu{1}.x;
            data_d_T{5}.x=data_Nu{1}.x;
            data_d_T{6}.x=data_Nu{1}.x;
            for data_ind=1:6
                [data_d_T_scaling(data_ind,1),data_d_T_scaling(data_ind,2)]=scaling(data_d_T{data_ind}.x(50:end),data_d_T{data_ind}.y(50:end));
            end
            data_d_T_sparse=sparse_data(data_d_T,4,5);
            
            plot_config.ylim_list=[1,0.0001,10000];
            plot_config.fontsize_legend=18;
            plot_config.xlim_list=[1,100,1000000];
            plot_config.ytick_list=[1,0.0001,0.001,0.01,0.1,1,10,100,1000];
            plot_config.loglog=[1,1];
            plot_config.legend_list={1,'$Nu$','$\textrm{max}_z[\partial_z \bar{T}_0(z)-1]$',...
                '$\textrm{max}_z[T_{\textrm{rms}}(z)]$','$1-\partial_z \bar{T}_0(z=0.5)$','$z_{\textrm{os}}$','$z_{p,T}$'};
            plot_config.user_color_style_marker_list={'k-','r--','b-.','ksquare','rx','bo','k-','r--','b-.'};    
            plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'Ra_overshoot.png'];
            plot_line(data_d_T_sparse,plot_config);
            
        case 'hewitt_3D_Ra'
            clear data_Nu
            %Read the Raleigh number and order my data by Rayleigh number
            for content_ind=1:length(dedalus_post_my)
                data_Nu{1}.x(content_ind)=dedalus_post_my{1,content_ind}.Ra_T;
            end
            [data_Nu{1}.x,ind]=sort(data_Nu{1}.x);
            dedalus_post_my=dedalus_post_my(1,ind);
            %---------read the Nusselt number.
            for content_ind=1:length(dedalus_post_my)
                data_Nu{1}.y(content_ind)=dedalus_post_my{1,content_ind}.Nu(1);
            end
            
            %%read the DNS data for comparison
            porous_otero_bound=dedalus_post_my{1}.get_porous_otero_bound();
            porous_hewitt_3D=dedalus_post_my{1}.get_porous_hewitt_3D();
            porous_pirozzoli_3D=dedalus_post_my{1}.get_porous_pirozzoli_3D();
            sparse_ind=4;
            data_Nu{2}.x=porous_hewitt_3D.Ra_Nu(1:sparse_ind:length(porous_hewitt_3D.Ra_Nu),1);
            data_Nu{2}.y=porous_hewitt_3D.Ra_Nu(1:sparse_ind:length(porous_hewitt_3D.Ra_Nu),2);
            data_Nu{3}.x=porous_pirozzoli_3D.Ra_Nu_DNS(1:length(porous_pirozzoli_3D.Ra_Nu_DNS),1);
            data_Nu{3}.y=porous_pirozzoli_3D.Ra_Nu_DNS(1:length(porous_pirozzoli_3D.Ra_Nu_DNS),2);
            sparse_ind=1;
            data_Nu{5}.x=porous_hewitt_3D.Ra_Nu(:,1);
            data_Nu{5}.y=porous_hewitt_3D.Ra_Nu(:,2);           
            data_Nu{6}.x=porous_pirozzoli_3D.Ra_Nu_DNS(:,1);
            data_Nu{6}.y=porous_pirozzoli_3D.Ra_Nu_DNS(:,2);
            data_Nu{4}.x=porous_otero_bound.Ra_Nu_bound(:,1);
            data_Nu{4}.y=porous_otero_bound.Ra_Nu_bound(:,2);
            
            %setup the plot config
            plot_config.linwidth=3;
            plot_config.label_list={1,'$Ra$','$Nu$'};
            plot_config.loglog=[1,1];
            plot_config.Markerindex=3;
            plot_config.xlim_list=0; plot_config.ylim_list=0;
            plot_config.print_size=[1,1000,1000];
            plot_config.user_color_style_marker_list={'k-','r^','rsquare','b--','r-','r-'}; 
            plot_config.legend_list={1,'Single mode','DNS (Hewitt \it{et al.} \rm 2014)','DNS (Pirozzoli \it{et al.} \rm 2021)','Upper bound (Otero \it{et al.} \rm 2004)'};
            plot_config.fontsize_legend=20;
            plot_config.fontsize=32;
            plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'Ra_Nu.png'];
            plot_config.linewidth=3;
            plot_config.ylim_list=[1,1,2500];
            plot_config.xtick_list=[1,1,10,100,1000,10000,100000,1000000];
            plot_config.ytick_list=[1,1,10,100,1000];
            plot_config.xlim_list=[1,10,1000000];
            plot_line(data_Nu,plot_config);
            [data_Nu_scaling.eta,data_Nu_scaling.c0]=scaling(data_Nu{1}.x(50:end),data_Nu{1}.y(50:end));
            %---------read RMS value
            for content_ind=1:length(dedalus_post_my)
                data_rms{1}.y(content_ind)=dedalus_post_my{1,content_ind}.T_rms_mid;
                data_rms{2}.y(content_ind)=dedalus_post_my{1,content_ind}.w_rms_mid;
                data_rms{3}.y(content_ind)=dedalus_post_my{1,content_ind}.u_rms_mid;
            end
            data_rms{1}.x=data_Nu{1}.x;
            data_rms{2}.x=data_Nu{1}.x;
            data_rms{3}.x=data_Nu{1}.x;
            data_rms{4}.x=porous_hewitt_3D.Ra_T_rms(:,1);
            data_rms{4}.y=porous_hewitt_3D.Ra_T_rms(:,2);
            data_rms{5}.x=porous_hewitt_3D.Ra_w_rms(:,1);
            data_rms{5}.y=porous_hewitt_3D.Ra_w_rms(:,2);
            data_rms{6}.x=porous_hewitt_3D.Ra_u_rms(:,1);
            data_rms{6}.y=porous_hewitt_3D.Ra_u_rms(:,2);
            %plot config of rms value
            plot_config.label_list={1,'$Ra$',[]};
            plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'Ra_rms.png'];
            plot_config.Markerindex=3;
            plot_config.user_color_style_marker_list={'r-','b--','k-.','rsquare','bo','k^'};    
            plot_config.loglog=[1,0];
            plot_config.linewidth=3;
            plot_config.legend_list={1,'${T_{\textrm{rms}}}(z=0.5)$','${w_{\textrm{rms}}}(z=0.5)$','${u_{\textrm{rms}}}(z=0.5)$',...
                '${T_{\textrm{rms}}}(z=0.5)$ (DNS)','${w_{\textrm{rms}}}(z=0.5)$ (DNS)','${u_{\textrm{rms}}}(z=0.5)$ (DNS)'};
            plot_config.fontsize_legend=18;
            plot_config.fontsize=32;
            plot_config.ylim_list=[1,0,0.15];
            plot_config.xlim_list=[1,100,1000000];
            plot_config.ytick_list=[1,0,0.05,0.1,0.15];
            plot_line(data_rms,plot_config);
            
            
            %---------read the mean gradient.  
            for content_ind=1:length(dedalus_post_my)
                d_T_total=dedalus_post_my{1,content_ind}.d_T_0+dedalus_post_my{1,content_ind}.dy_T_mean;
                data_d_T{1}.y(content_ind)=dedalus_post_my{1,content_ind}.Nu(1);
                data_d_T{2}.y(content_ind)=max(d_T_total);
                data_d_T{3}.y(content_ind)=dedalus_post_my{1,content_ind}.T_rms_max;
                data_d_T{4}.y(content_ind)=-dedalus_post_my{1,content_ind}.d_T_0_mid;
                data_d_T{5}.y(content_ind)=dedalus_post_my{1,content_ind}.z_T_BL;
                data_d_T{6}.y(content_ind)=dedalus_post_my{1,content_ind}.z_T_rms_max;
                if data_d_T{6}.y(content_ind)>0.5
                    data_d_T{6}.y(content_ind)=1-data_d_T{6}.y(content_ind);
                end
            end
            negative_max_d_T_0=find(data_d_T{2}.y<0);
            data_d_T{5}.y(negative_max_d_T_0)=NaN;
            
            data_d_T{1}.x=data_Nu{1}.x;
            data_d_T{2}.x=data_Nu{1}.x;
            data_d_T{3}.x=data_Nu{1}.x;
            data_d_T{4}.x=data_Nu{1}.x;
            data_d_T{5}.x=data_Nu{1}.x;
            data_d_T{6}.x=data_Nu{1}.x;
            for data_ind=1:6
                [data_d_T_scaling(data_ind,1),data_d_T_scaling(data_ind,2)]=scaling(data_d_T{data_ind}.x(50:end),data_d_T{data_ind}.y(50:end));
            end
            data_d_T_sparse=sparse_data(data_d_T,4,5);
            
            plot_config.ylim_list=[1,0.0001,10000];
            plot_config.fontsize_legend=18;
            plot_config.ytick_list=[1,0.0001,0.001,0.01,0.1,1,10,100,1000];
            plot_config.loglog=[1,1];
            plot_config.legend_list={1,'$Nu$','${\textrm{max}}_z[\partial_z \bar{T}_0(z)-1]$',...
                '${\textrm{max}}_z[T_{\textrm{rms}}(z)]$','$1-\partial_z \bar{T}_0(z=0.5)$','$z_{\textrm{os}}$','$z_{p,T}$'};
            plot_config.user_color_style_marker_list={'k-','r--','b-.','ksquare','rx','bo','k-','r--','b-.'};    
            plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'Ra_overshoot.png'];
            plot_line(data_d_T_sparse,plot_config);
        
        case 'HB_porous_Nu_kx_Ra'
            for Ra_ind=1:size(dedalus_post_my,1)
                for kx_ind=1:size(dedalus_post_my,2)
                    data{Ra_ind}.y(kx_ind)=dedalus_post_my{Ra_ind,kx_ind}.Nu(1);
                    data{Ra_ind}.x(kx_ind)=dedalus_post_my{Ra_ind,kx_ind}.kx;
                    [data{Ra_ind}.x,ind]=sort(data{Ra_ind}.x);
                    data{Ra_ind}.y=data{Ra_ind}.y(ind);
                    dedalus_post_my(Ra_ind,ind)=dedalus_post_my(Ra_ind,ind);
                end
                Ra_list(Ra_ind)=dedalus_post_my{Ra_ind,1}.Ra_T;
            end
            [Ra_list,ind]=sort(Ra_list);
            data=data(ind);
            for Ra_ind=1:length(Ra_list)
               %ADD the wavenumber that maximize the 
               [val,kx_max_ind]=max(data{Ra_ind}.y);
               data{11}.x(Ra_ind)=data{Ra_ind}.x(kx_max_ind);
               data{11}.y(Ra_ind)=val;
               Ra=Ra_list(Ra_ind);
               data_rms_Ra{1}.x(Ra_ind)=Ra_list(Ra_ind);
               data_rms_Ra{1}.y(Ra_ind)=dedalus_post_my{Ra_ind,kx_ind}.T_rms_mid;
               data_rms_Ra{2}.x(Ra_ind)=Ra_list(Ra_ind);
               data_rms_Ra{2}.y(Ra_ind)=dedalus_post_my{Ra_ind,kx_ind}.w_rms_mid;
               data_rms_Ra{3}.x(Ra_ind)=Ra_list(Ra_ind);
               data_rms_Ra{3}.y(Ra_ind)=dedalus_post_my{Ra_ind,kx_ind}.u_rms_mid;
               %Add the wavenumber that I am using from Hewitt.. and
               %correspoding Nusselt number.
               data{12}.x(Ra_ind)=0.48*Ra^(0.4);
               [kx_val,kx_ind]=min(abs(data{12}.x(Ra_ind)-data{Ra_ind}.x));
               %data{12}.x(Ra_ind)=data{Ra_ind}.x(kx_ind);
               data{12}.y(Ra_ind)=data{Ra_ind}.y(kx_ind);
            end
            %Add the wavenumber corresponding to the Nusselt number in the
            %DNS, Ra=5000
            Ra_ind=1; Nu=37.249;
            data{13}.x=linspace(1,200,200);
            data{13}.y=ones(1,200)*Nu;
%             [min]=min(data{Ra_ind}.y-Nu));
            
            %Add the wavenumber corresponding to the Nusselt number in the
            %DNS, Ra=10000
            Ra_ind=3; Nu=72.198;
%             data{14}.x=linspace(1,200,200);
%             data{14}.y=ones(1,200)*Nu;
%             
            %Add the wavenumber corresponding to the Nusselt number in the
            %DNS, Ra=20000
            Ra_ind=7; Nu=142.457;
%             data{15}.x=linspace(1,200,200);
%             data{15}.y=ones(1,200)*Nu;
%             
            %Add the wavenumber corresponding to the Nusselt number in the
            %DNS, Ra=40000
            Ra_ind=10; Nu=281.086;
            data{14}.x=linspace(1,200,200);
            data{14}.y=ones(1,200)*Nu;
            data{15}.x=10;
            data{15}.y=37.25;
            plot_config.user_color_style_marker_list={'k-','k-','k-','k-','k-','k-','k-','k-','k-','k-','b*','ro','m--','m--','msquare','m--x'};
            plot_config.Markerindex=3;
            plot_config.xtick_list=[1,1,10,100];
            plot_config.ytick_list=[1,1,10,100,1000,10000];
            plot_config.label_list={1,'$k_x$','$Nu$'};
            plot_config.loglog=[1,1];
            plot_config.print_size=[1,1000,1000];
            plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'Nu_kx.png'];
            plot_config.fontsize=32;
            plot_config.linewidth=3;
            plot_line(data,plot_config);
           
        case 'HB_porous_kx'
            for kx_ind=1:length(dedalus_post_my)
                kx_list(kx_ind)=dedalus_post_my{kx_ind}.kx;
            end
            [kx_list,ind]=sort(kx_list);
            dedalus_post_my=dedalus_post_my(ind);
            
            for ind=1:8
                kx_plot(ind)=dedalus_post_my{ind*10}.kx;
                data{ind}.x=dedalus_post_my{ind*10}.T_0+1-dedalus_post_my{ind*5}.z_list;
                data{ind}.y=dedalus_post_my{ind*10}.z_list;
            end
            data=sparse_data(data,5,200);
            data{13}.x=NaN; data{13}.y=NaN;
            data{14}.x=NaN; data{14}.y=NaN;
            data{15}.x=NaN; data{15}.y=NaN;
            data{16}.x=NaN; data{16}.y=NaN;
            plot_config.legend_index=[1,2,3,4,13,14,15,16];
            plot_config.user_color_style_marker_list={'k-','b--','r-.','m:','ko','bsquare','r*','mx','k-','b--','r-.','m:','k-o','b--square','r-.*','m:x'};
            plot_config.legend_list={1,'$k_x=10$','$k_x=20$','$k_x=30$','$k_x=40$','$k_x=50$','$k_x=60$','$k_x=70$','$k_x=80$'};
            plot_config.label_list={1,'$\bar{T}_0+1-z$','$z$'};
            plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'T_0_kx.png'];
            plot_config.xlim_list=[1,0,1];
            plot_config.ylim_list=[1,0,1];
            plot_config.print_size=[1,500,1000];
            plot_config.fontsize_legend=18;
            plot_config.Markerindex=3;   
            plot_config.fontsize=32;
            plot_config.linewidth=3;
            plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'T_0_kx.png'];
            plot_line(data,plot_config);
            plot_config.legend_list={0};
            plot_config.ylim_list=[1,0,0.03];
            plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'T_0_kx_boundary.png'];
            plot_line(data,plot_config);
%             
%             
        case 'hewitt_2_layer_Omega'
            %Read the Raleigh number and order my data by Rayleigh number
            for content_ind=1:length(dedalus_post_my)
                Omega_list(content_ind)=dedalus_post_my{1,content_ind}.HB_porous_2_layer_Omega;
            end
            [Omega_list,ind]=sort(Omega_list);
            dedalus_post_my=dedalus_post_my(1,ind);
            %read the DNS data from Hewitt...
            porous_hewitt_2_layer=dedalus_post_my{1}.get_porous_hewitt_2_layer();

            %---------read Nusselt number and compare
            for content_ind=1:length(dedalus_post_my)
                data_Nu{1}.y(content_ind)=dedalus_post_my{1,content_ind}.Nu(1);
            end
            data_Nu{1}.x=Omega_list;
            data_Nu{2}.x=porous_hewitt_2_layer.Omega_Nu(:,1);
            data_Nu{2}.y=porous_hewitt_2_layer.Omega_Nu(:,2);
            
            %setup the plot config
            plot_config.linwidth=3;
            plot_config.label_list={1,'$\Omega$','$Nu$'};
            plot_config.loglog=[1,1];
            plot_config.Markerindex=3;
            plot_config.xlim_list=0; plot_config.ylim_list=0;
            plot_config.print_size=[1,1000,1000];
            plot_config.user_color_style_marker_list={'k-','ro','b--','r-'}; 
            plot_config.legend_list={1,'Single mode','DNS (Hewitt \it{et al.} \rm 2014)'};
            plot_config.fontsize_legend=20;
            plot_config.fontsize=32;
            plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'Ra_Nu.png'];
            plot_config.linewidth=3;
            plot_config.ylim_list=[1,10,50];
            plot_config.xtick_list=[1,0.001,0.01,0.1,1,10,100];
            plot_config.ytick_list=[1,10,20,30,40];
            plot_config.xlim_list=[1,0.001,100];
            plot_config.loglog=[1,0];
            plot_line(data_Nu,plot_config);
            
            
            %---------read RMS value and plot their trend over Omega
            for content_ind=1:length(dedalus_post_my)
                data_rms{1}.y(content_ind)=dedalus_post_my{1,content_ind}.T_rms_mid;
                data_rms{2}.y(content_ind)=dedalus_post_my{1,content_ind}.w_rms_mid;
                data_rms{3}.y(content_ind)=dedalus_post_my{1,content_ind}.u_rms_mid;
                data_rms{1}.x(content_ind)=dedalus_post_my{1,content_ind}.HB_porous_2_layer_Omega;
            end
            data_rms{2}.x=data_rms{1}.x;
            data_rms{3}.x=data_rms{1}.x;
            data_rms{4}.x=porous_hewitt_2_layer.Omega_T_rms(:,1);
            data_rms{4}.y=porous_hewitt_2_layer.Omega_T_rms(:,2);
            data_rms{5}.x=porous_hewitt_2_layer.Omega_w_rms(:,1);
            data_rms{5}.y=porous_hewitt_2_layer.Omega_w_rms(:,2);
            data_rms{6}.x=porous_hewitt_2_layer.Omega_u_rms(:,1);
            data_rms{6}.y=porous_hewitt_2_layer.Omega_u_rms(:,2);
            %plot config of rms value
            plot_config.label_list={1,'$\Omega$',[]};
            plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'Omega_rms.png'];
            plot_config.Markerindex=3;
            plot_config.user_color_style_marker_list={'r-','b--','k-.','rsquare','bo','k^'};    
            plot_config.loglog=[1,0];
            plot_config.linewidth=3;
            plot_config.legend_list={1,'$T_{\textrm{rms}}(z=0.5)$','$w_{\textrm{rms}}(z=0.5)$','$u_{\textrm{rms}}(z=0.5)$',...
                '${T_{\textrm{rms}}}(z=0.5)$ (DNS)','${w_{\textrm{rms}}}(z=0.5)$ (DNS)','$u_{\textrm{rms}}(z=0.5)$ (DNS)'};
            plot_config.fontsize_legend=24;
            plot_config.fontsize=32;
            plot_config.ylim_list=[1,0,0.25];
            plot_config.xlim_list=[1,0.001,100];
            plot_config.xtick_list=[1,0.001,0.01,0.1,1,10,100];
            plot_config.ytick_list=[1,0,0.05,0.1,0.15,0.2,0.25];
            plot_line(data_rms,plot_config);
            
            %---------plot the profile over z for T_0 and rms.. for
            %comparison
            %get the index that corresponds to the Omega=0.04,0.25,1.28,10
            [~,Omega_ind(1)]=min(abs(Omega_list-0.04));
            [~,Omega_ind(2)]=min(abs(Omega_list-0.25));
            [~,Omega_ind(3)]=min(abs(Omega_list-1.28));
            [~,Omega_ind(4)]=min(abs(Omega_list-10));
            
            for ind=1:length(Omega_ind)
                background_T=1-dedalus_post_my{Omega_ind(ind)}.z_list;
                data_T{ind}.y=dedalus_post_my{Omega_ind(ind)}.z_list;
                data_T{ind}.x=dedalus_post_my{Omega_ind(ind)}.T_0+background_T;
            end
            data_T{5}.x=porous_hewitt_2_layer.z_T_0_Omega_0p04(:,2);
            data_T{5}.y=porous_hewitt_2_layer.z_T_0_Omega_0p04(:,1);
            data_T{6}.x=porous_hewitt_2_layer.z_T_0_Omega_0p25(:,2);
            data_T{6}.y=porous_hewitt_2_layer.z_T_0_Omega_0p25(:,1);
            data_T{7}.x=porous_hewitt_2_layer.z_T_0_Omega_1p28(:,2);
            data_T{7}.y=porous_hewitt_2_layer.z_T_0_Omega_1p28(:,1);
            data_T{8}.x=porous_hewitt_2_layer.z_T_0_Omega_10(:,2);
            data_T{8}.y=porous_hewitt_2_layer.z_T_0_Omega_10(:,1);
            data_T=sparse_data(data_T,5,5);
            
            plot_config.label_list={1,'$\bar{T}_0+1-z$','$z$'};
            plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'T_0.png'];
            plot_config.legend_list={1,'$\Omega=0.04$','$\Omega=0.25$', '$\Omega=1.28$','$\Omega=10$'...
                    '$\Omega=0.04$ (DNS)', '$\Omega=0.25$ (DNS)', '$\Omega=1.28$ (DNS)','$\Omega=10$ (DNS)'};
            plot_config.fontsize_legend=22;
            plot_config.fontsize=32;
            plot_config.linwidth=3;
            plot_config.Markerindex=3;
            plot_config.user_color_style_marker_list={'r-','k--','b-.','m:','rsquare','k^','bo','mx','r-','k--','b-.','m:'};
            plot_config.print_size=[1,900,1000];
            plot_config.xlim_list=[1,0,1];
            plot_config.ylim_list=[1,0,1];
            plot_config.loglog=[0,0];
            plot_config.xtick_list=[1,0,0.2,0.4,0.6,0.8,1];
            plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
            plot_line(data_T,plot_config);
           
            %plot of T_rms comparison...
            rms_list={'T','w','u'};
            plot_config.linwidth=3;
            plot_config.fontsize=32;
            plot_config.legend_list={0};
            for rms_ind=1:length(rms_list)
                clear data_T
                rms_name=rms_list{rms_ind};
                if rms_name=='u'
                    suffix='_tilde';
                else 
                    suffix='_hat';
                end
                data_T{1}.x=abs(dedalus_post_my{Omega_ind(1)}.([rms_name,suffix]))*sqrt(2);
                data_T{1}.y=dedalus_post_my{Omega_ind(1)}.z_list;
                data_T{2}.x=abs(dedalus_post_my{Omega_ind(2)}.([rms_name,suffix]))*sqrt(2);
                data_T{2}.y=dedalus_post_my{Omega_ind(2)}.z_list;
                data_T{3}.x=abs(dedalus_post_my{Omega_ind(3)}.([rms_name,suffix]))*sqrt(2);
                data_T{3}.y=dedalus_post_my{Omega_ind(3)}.z_list;
                data_T{4}.x=abs(dedalus_post_my{Omega_ind(4)}.([rms_name,suffix]))*sqrt(2);
                data_T{4}.y=dedalus_post_my{Omega_ind(4)}.z_list;

                porous_hewitt_3D=dedalus_post_my{1}.get_porous_hewitt_3D();
                data_T{5}.x=porous_hewitt_2_layer.(['z_',rms_name,'_rms_Omega_0p04'])(:,2);
                data_T{5}.y=porous_hewitt_2_layer.(['z_',rms_name,'_rms_Omega_0p04'])(:,1);
                data_T{6}.x=porous_hewitt_2_layer.(['z_',rms_name,'_rms_Omega_0p25'])(:,2);
                data_T{6}.y=porous_hewitt_2_layer.(['z_',rms_name,'_rms_Omega_0p25'])(:,1);
                data_T{7}.x=porous_hewitt_2_layer.(['z_',rms_name,'_rms_Omega_1p28'])(:,2);
                data_T{7}.y=porous_hewitt_2_layer.(['z_',rms_name,'_rms_Omega_1p28'])(:,1);
                data_T{8}.x=porous_hewitt_2_layer.(['z_',rms_name,'_rms_Omega_10'])(:,2);
                data_T{8}.y=porous_hewitt_2_layer.(['z_',rms_name,'_rms_Omega_10'])(:,1);
                %data_T=sparse_data(data_T,5,5);

                %                 sparse_ind=1;
%                 data_T{4}.x=porous_hewitt_3D.(['Ra_4000_z_',rms_name,'_rms'])(1:sparse_ind:length(porous_hewitt_3D.(['Ra_4000_z_',rms_name,'_rms'])),2);
%                 data_T{4}.y=porous_hewitt_3D.(['Ra_4000_z_',rms_name,'_rms'])(1:sparse_ind:length(porous_hewitt_3D.(['Ra_4000_z_',rms_name,'_rms'])),1);
%                 data_T{5}.x=porous_hewitt_3D.(['Ra_8000_z_',rms_name,'_rms'])(1:sparse_ind:length(porous_hewitt_3D.(['Ra_8000_z_',rms_name,'_rms'])),2);
%                 data_T{5}.y=porous_hewitt_3D.(['Ra_8000_z_',rms_name,'_rms'])(1:sparse_ind:length(porous_hewitt_3D.(['Ra_8000_z_',rms_name,'_rms'])),1);
%                 data_T{6}.x=porous_hewitt_3D.(['Ra_16000_z_',rms_name,'_rms'])(1:sparse_ind:length(porous_hewitt_3D.(['Ra_16000_z_',rms_name,'_rms'])),2);
%                 data_T{6}.y=porous_hewitt_3D.(['Ra_16000_z_',rms_name,'_rms'])(1:sparse_ind:length(porous_hewitt_3D.(['Ra_16000_z_',rms_name,'_rms'])),1);
% 
%                 sparse_ind=5;
%                 data_T{7}.x=porous_hewitt_3D.(['Ra_4000_z_',rms_name,'_rms'])(1:sparse_ind:length(porous_hewitt_3D.(['Ra_4000_z_',rms_name,'_rms'])),2);
%                 data_T{7}.y=porous_hewitt_3D.(['Ra_4000_z_',rms_name,'_rms'])(1:sparse_ind:length(porous_hewitt_3D.(['Ra_4000_z_',rms_name,'_rms'])),1);
%                 data_T{8}.x=porous_hewitt_3D.(['Ra_8000_z_',rms_name,'_rms'])(1:sparse_ind:length(porous_hewitt_3D.(['Ra_8000_z_',rms_name,'_rms'])),2);
%                 data_T{8}.y=porous_hewitt_3D.(['Ra_8000_z_',rms_name,'_rms'])(1:sparse_ind:length(porous_hewitt_3D.(['Ra_8000_z_',rms_name,'_rms'])),1);
%                 data_T{9}.x=porous_hewitt_3D.(['Ra_16000_z_',rms_name,'_rms'])(1:sparse_ind:length(porous_hewitt_3D.(['Ra_16000_z_',rms_name,'_rms'])),2);
%                 data_T{9}.y=porous_hewitt_3D.(['Ra_16000_z_',rms_name,'_rms'])(1:sparse_ind:length(porous_hewitt_3D.(['Ra_16000_z_',rms_name,'_rms'])),1);

                plot_config.label_list={1,['$',rms_name,'_{\textrm{rms}}(z)$'],'$z$'};
                plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,rms_name,'_rms.png'];
                plot_config.Markerindex=3;
                plot_config.user_color_style_marker_list={'r-','k--','b-.','m:','rsquare','k^','bo','mx','r-','k--','b-.','m:'};
                plot_config.print_size=[1,600,1000];
                switch rms_name
                    case 'T'
                        plot_config.xlim_list=[1,0,0.3];
                    case 'u'
                        plot_config.xlim_list=[1,0,0.2];
                    case 'w'
                        plot_config.xlim_list=[1,0,0.1];
                end
                plot_config.xtick_list=[1,0,0.1,0.2,0.3];
                plot_line(data_T,plot_config);
            end
            
        case 'HB_porous_thermal_BC'
            plot_config.label_list={1,'$\bar{T}_0+1-z$','$z$'};
            plot_config.Markerindex=3;
            plot_config.print_size=[1,600,1000];
            marker_num=5;
            marker_ind=round(linspace(1,length(data_T{1}.x),marker_num));
            marker_ind=marker_ind(2:end-1);
            for i=2:5
                data_T{i+4}.x=data_T{i}.x(marker_ind);
                data_T{i+4}.y=data_T{i}.y(marker_ind);
            end
            plot_config.user_color_style_marker_list={'k-','b--','b--','r--','r--','bsquare','bo','r^','rv'};
            plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'T_0.png'];
            plot_config.print=1;
            plot_line(data_T,plot_config);
        case {'HB_porous_diffusive_BC','HB_porous_finger_BC'}
            plot_config.Markerindex=3;
            plot_config.print_size=[1,600,1000];
            marker_num=5;
            marker_ind=round(linspace(1,length(data_T{1}.x),marker_num));
            marker_ind=marker_ind(2:end-1);
            for i=2:7
                data_T{i+6}.x=data_T{i}.x(marker_ind);
                data_T{i+6}.y=data_T{i}.y(marker_ind);
                data_S{i+6}.x=data_S{i}.x(marker_ind);
                data_S{i+6}.y=data_S{i}.y(marker_ind);
            end
            plot_config.user_color_style_marker_list={'k-','b--','b--','r--','r--','m--','m--','bsquare','bo','r^','rv','mx','m+'};
            plot_config.print=1;
            plot_config.label_list={1,'$\bar{T}_0+1-z$','$z$'};
            plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'T_0.png'];
            plot_line(data_T,plot_config);
            plot_config.label_list={1,'$\bar{S}_0+1-z$','$z$'};
            plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'S_0.png'];
            plot_line(data_S,plot_config);


        case 'HB_porous_finger_BC'
        case 'HB_benard_kx'

            Nu_kx=dedalus_post_my{1}.get_Nu_kx_Toomre();
            [data_Nu{1}.x,ind]=sort(data_Nu{1}.x);
            data_Nu{1}.y=data_Nu{1}.y(ind);
            data_Nu{2}.x=Nu_kx(:,1);
            data_Nu{2}.y=Nu_kx(:,2);
            plot_config.Markerindex=3;
            plot_config.user_color_style_marker_list={'b-','k*'};
            plot_config.label_list={1,'$k_x$','Nu'};
            plot_config.print_size=[1,1000,900];
            plot_config.name=['C:\Figure\DDC_LST\',group_name,'Nu_kx_Toomre_1977.png'];
            plot_line(data_Nu,plot_config);


            clear plot_config;
            T_mean_var=1-dedalus_post_my{1}.z_list;
            T_mean_sign='+1-z';
            plot_config.legend_list{1}=1;
            y_ind=round(linspace(1,length(dedalus_post_my{1}.z_list),5));
            for kx_ind=1:6
                data_T0{kx_ind}.x=dedalus_post_my{2+5*(kx_ind-1)}.T_0+T_mean_var;
                data_T0{kx_ind}.y=dedalus_post_my{2+5*(kx_ind-1)}.z_list;
                plot_config.legend_list{kx_ind+1}=['$k_x=',num2str(round(dedalus_post_my{2+5*(kx_ind-1)}.kx)),'$'];
                data_T0{kx_ind+6}.x=data_T0{kx_ind}.x(y_ind);
                data_T0{kx_ind+6}.y=data_T0{kx_ind}.y(y_ind);
            end
            plot_config.ylim_list=[1,0,1];
            plot_config.xlim_list=[1,0,1];
            plot_config.print_size=[1,900,900];
            plot_config.loglog=[0,0];
            plot_config.label_list={1,['$\bar{T}_0',T_mean_sign,'$'], '$z$'};
            plot_config.Markerindex=3;
            plot_config.user_color_style_marker_list={'k-','k-','b--','b--','r-.','r-.'...
                ,'ksquare','ko','b^','bv','rx','r+'};
            plot_config.fontsize_legend=20;
            plot_config.name=['C:\Figure\DDC_LST\',group_name,'T_0_kx.png'];
            plot_line(data_T0,plot_config);


        case 'HB_benard_Ra'

            Nu_Ra=dedalus_post_my{1}.get_Nu_Ra_Toomre();
            [data_Nu{1}.x,ind]=sort(data_Nu{1}.x);
            data_Nu{1}.y=data_Nu{1}.y(ind);
            data_Nu{2}.x=Nu_Ra(:,1);
            data_Nu{2}.y=Nu_Ra(:,2)+1;
            plot_config.ylim_list=[1,1,100];
            plot_config.ytick_list=[1,1,10,100];
            plot_config.xtick_list=[1,10^5,10^6,10^7,10^8,10^9,10^10];
            plot_config.Markerindex=3;
            plot_config.loglog=[1,1];
            plot_config.user_color_style_marker_list={'b-','k*'};
            plot_config.label_list={1,'$Ra$','Nu'};
            plot_config.print_size=[1,1000,900];
            plot_config.name=['C:\Figure\DDC_LST\',group_name,'Nu_Ra_Toomre_1977.png'];
            plot_line(data_Nu,plot_config);

            clear plot_config;
            T_mean_var=1-dedalus_post_my{1}.z_list;
            T_mean_sign='+1-z';
            plot_config.legend_list{1}=1;
            y_ind=round(linspace(1,length(dedalus_post_my{1}.z_list),5));
            for Ra_ind=1:6
                data_T0{Ra_ind}.x=dedalus_post_my{1+5*(Ra_ind-1)}.T_0+T_mean_var;
                data_T0{Ra_ind}.y=dedalus_post_my{1+5*(Ra_ind-1)}.z_list;
                plot_config.legend_list{Ra_ind+1}=['$Ra_T=',num2str(round(dedalus_post_my{1+5*(Ra_ind-1)}.Ra_T)),'$'];
                data_T0{Ra_ind+6}.x=data_T0{Ra_ind}.x(y_ind);
                data_T0{Ra_ind+6}.y=data_T0{Ra_ind}.y(y_ind);
            end
            plot_config.ylim_list=[1,0,1];
            plot_config.xlim_list=[1,0,1];
            plot_config.print_size=[1,900,900];
            plot_config.loglog=[0,0];
            plot_config.label_list={1,['$\bar{T}_0',T_mean_sign,'$'], '$z$'};
            plot_config.Markerindex=3;
            plot_config.user_color_style_marker_list={'k-','k-','b--','b--','r-.','r-.'...
                ,'ksquare','ko','b^','bv','rx','r+'};
            plot_config.fontsize_legend=20;
            plot_config.name=['C:\Figure\DDC_LST\',group_name,'T_0_Ra.png'];
            plot_line(data_T0,plot_config);
            
        case 'HB_benard_yang'
            z_list=dedalus_post_my{1,1}.z_list;
            for slurm_ind=1:size(dedalus_post_my,1)
                for content_ind=1:size(dedalus_post_my,2)
                    Ra_S2T(slurm_ind,content_ind)=dedalus_post_my{slurm_ind,content_ind}.Ra_S2T;
                    S_0(:,content_ind)=dedalus_post_my{slurm_ind,content_ind}.S_0;
                    T_0(:,content_ind)=dedalus_post_my{slurm_ind,content_ind}.T_0;
                end
                [Ra_S2T_tmp,ind]=sort(Ra_S2T(slurm_ind,:));
                Ra_S2T(slurm_ind,:)=Ra_S2T_tmp;  
                S_0=S_0(:,ind);
                T_0=T_0(:,ind);
                
                content_len=size(dedalus_post_my,2);
                for content_ind=1:size(dedalus_post_my,2)
                    data{content_ind}.x=S_0(:,content_ind)+z_list+content_ind-1;
                    data{content_ind}.y=z_list;
                    data{content_ind+content_len}.x=T_0(:,content_ind)+z_list+content_ind-1;
                    data{content_ind+content_len}.y=z_list;
                end
                plot_config.label_list={1,'','$z$'};
                plot_config.Markerindex=3;
                plot_config.user_color_style_marker_list={'k-','k-','k-','k-','k-','k-','k-','b--','b--','b--','b--','b--','b--','b--'};
                switch slurm_ind
                    case 1
                         plot_config.name=['C:\Figure\DDC_LST\',group_name,'_Ra_1e5.png'];
                    case 2
                         plot_config.name=['C:\Figure\DDC_LST\',group_name,'_Ra_1e6.png'];
                end
                plot_config.xlim_list=[1,-0.1,7.1];
                plot_config.ylim_list=[1,-0.1,1.1];
                plot_config.xtick_list=[1,0,1,2,3,4,5,6,7];
                plot_config.ytick_list=[1,0,0.5,1];
                plot_line(data,plot_config);
            end
            
        case 'HB_benard_salt_finger_kx'
            for slurm_ind=1:size(dedalus_post_my,1)
                for content_ind=1:size(dedalus_post_my,2)
                    if ~isempty(dedalus_post_my{slurm_ind,content_ind})
                        data{slurm_ind}.x(content_ind)=dedalus_post_my{slurm_ind,content_ind}.kx;
                        data{slurm_ind}.y(content_ind)=abs(dedalus_post_my{slurm_ind,content_ind}.Nu_S(1));

                        eig_inf_ind=find(real(dedalus_post_my{slurm_ind,content_ind}.eigenvalues)>10^6);
                        dedalus_post_my{slurm_ind,content_ind}.eigenvalues(eig_inf_ind)=[];
                        lambda_max{slurm_ind}.y(content_ind)=dedalus_post_my{slurm_ind,content_ind}.eigenvalues(1);
                    end
                end
            end
            sol_num=7;
            for sol_ind=1:sol_num
                %%combine the data that increase R_rho_S2T and decrease
                %%R_rho_S2T
                data_all{sol_ind}.x=[data{sol_ind}.x,data{sol_ind+sol_num}.x];
                data_all{sol_ind}.y=[data{sol_ind}.y,data{sol_ind+sol_num}.y];
                lambda_max_all{sol_ind}.y=[lambda_max{sol_ind}.y,lambda_max{sol_ind+sol_num}.y];
                
                %sort data
                [~,ind]=sort(data_all{sol_ind}.x);
%                 ind=ind(2:end);
                data_all{sol_ind}.x=data_all{sol_ind}.x(ind);
                data_all{sol_ind}.y=data_all{sol_ind}.y(ind);
                lambda_max_all{sol_ind}.x=data_all{sol_ind}.x;
                lambda_max_all{sol_ind}.y=lambda_max_all{sol_ind}.y(ind);
                stable_ind=find(real(lambda_max_all{sol_ind}.y)<=0);
                if isempty(stable_ind)
                    data_all_stable{sol_ind}.x=NaN;
                    data_all_stable{sol_ind}.y=NaN;
                    data_all_stable_2{sol_ind}.x=NaN;
                    data_all_stable_2{sol_ind}.y=NaN;
                elseif isempty(find(diff(stable_ind)>1))
                    data_all_stable{sol_ind}.x=data_all{sol_ind}.x(stable_ind);
                    data_all_stable{sol_ind}.y=data_all{sol_ind}.y(stable_ind);
                    data_all_stable_2{sol_ind}.x=NaN; data_all_stable_2{sol_ind}.y=NaN;
                else
                    stable_2_ind=find(diff(stable_ind)>1);
                    stable_2_ind=stable_2_ind(1);
                    data_all_stable{sol_ind}.x=data_all{sol_ind}.x(stable_ind(1:stable_2_ind));
                    data_all_stable{sol_ind}.y=data_all{sol_ind}.y(stable_ind(1:stable_2_ind));
                    data_all_stable_2{sol_ind}.x=data_all{sol_ind}.x(stable_ind(1+stable_2_ind:end));
                    data_all_stable_2{sol_ind}.y=data_all{sol_ind}.y(stable_ind(1+stable_2_ind:end));
                end
                
                unstable_ind=find(real(lambda_max_all{sol_ind}.y)>0);
                if isempty(unstable_ind)
                    data_all_unstable{sol_ind}.x=NaN;
                    data_all_unstable{sol_ind}.y=NaN;
                    data_all_unstable_2{sol_ind}.x=NaN; data_all_unstable_2{sol_ind}.y=NaN;
                elseif isempty(find(diff(unstable_ind)>1, 1))
                    data_all_unstable{sol_ind}.x=data_all{sol_ind}.x(unstable_ind);
                    data_all_unstable{sol_ind}.y=data_all{sol_ind}.y(unstable_ind);
                    data_all_unstable_2{sol_ind}.x=NaN; data_all_unstable_2{sol_ind}.y=NaN;
                else
                    unstable_2_ind=find(diff(unstable_ind)>1);
                    unstable_2_ind=unstable_2_ind(1);
                    data_all_unstable{sol_ind}.x=data_all{sol_ind}.x(unstable_ind(1:unstable_2_ind));
                    data_all_unstable{sol_ind}.y=data_all{sol_ind}.y(unstable_ind(1:unstable_2_ind));
                    data_all_unstable_2{sol_ind}.x=data_all{sol_ind}.x(unstable_ind(1+unstable_2_ind:end));
                    data_all_unstable_2{sol_ind}.y=data_all{sol_ind}.y(unstable_ind(1+unstable_2_ind:end));
                end
%                 if all(real(lambda_max{sol_ind}.y)>0)
%                     %plot_config.user_color_style_marker_list{sol_ind}=[color_list{sol_ind},'--'];
%                 elseif all(real(lambda_max{sol_ind}.y)<0)
%                     %plot_config.user_color_style_marker_list{sol_ind}=[color_list{sol_ind},'-'];
%                 else 
%                     error('There is a transition of stability in parameter regime');
%                 end
            end
            %Nu_salt_finger_yang=dedalus_post_my{1,1}.get_Nu_salt_finger_yang();
            plot_config.legend_list={1,'S2','A2','S2L','A3','S1','A1','S3'};
            order=[5,1,7,3,6,2,4];%This will give S1, S2, S3, S2L, A1, A2, A3 order...
            plot_config.legend_list(2:end)=plot_config.legend_list(order+1);
            data_kx=[data_all_stable(order),data_all_stable_2(order),data_all_unstable(order),data_all_unstable_2(order)];
            plot_config.user_color_style_marker_list={'-','-','-','-','-','-','-',...
                                                    '-','-','-','-','-','-','-',...                                        
                                                    ':',':',':',':',':',':',':',...
                                                    ':',':',':',':',':',':',':'};
            data_len=length(data_kx);
            %data_R_rho_T2S{data_len+1}.x=Nu_salt_finger_yang.Ra_T_1e5_Nu_S(:,1);
            %data_R_rho_T2S{data_len+1}.y=Nu_salt_finger_yang.Ra_T_1e5_Nu_S(:,2);
            %plot_config.user_color_style_marker_list{data_len+1}='k*';
            plot_config.Markerindex=3;
            plot_config.label_list={1,'$k_x$','$Sh$'};
           
            plot_config.linewidth=5;
            plot_config.RGB={1,[0,0,0],[1,0,0],[0,0,1],[0,1,1],[1,0,1],...
                                [0.8500 0.3250 0.0980],[0,1,0],...
                                [0,0,0],[1,0,0],[0,0,1],[0,1,1],[1,0,1],...
                                [0.8500 0.3250 0.0980],[0,1,0],...
                                [0,0,0],[1,0,0],[0,0,1],[0,1,1],[1,0,1],...
                                [0.8500 0.3250 0.0980],[0,1,0],...
                                [0,0,0],[1,0,0],[0,0,1],[0,1,1],[1,0,1],...
                                [0.8500 0.3250 0.0980],[0,1,0],...
                                [0,0,0]};
            plot_config.name=['C:\Figure\DDC_LST\',group_name,'_Sh_kx.png'];
            plot_config.loglog=[0,0];
            plot_config.ylim_list=[1,1,25];
            plot_config.ytick_list=[1,1,5,10,15,20,25];
            plot_config.xlim_list=[0,0,60];
%             plot_config.xtick_list=[0,1,10,60];
            plot_line(data_kx,plot_config);
            
            
%             for content_ind=1:size(dedalus_post_my,2)
%                 Nu(content_ind)=dedalus_post_my{1,content_ind}.Nu(1);
%                 Nu_S(content_ind)=dedalus_post_my{1,content_ind}.Nu_S(1);
%                 kx(content_ind)=dedalus_post_my{1,content_ind}.kx;
%             end

        case 'HB_benard_salt_finger_profile'
            legend_list={'S2','A2','S2L','A3','S1','A1','S3'};
            order=[5,1,7,3,6,2,4];%This will give S1, S2, S3, S2L, A1, A2, A3 order...
            
            legend_list=legend_list(order);
            RGB={[0,0,0],[1,0,0],[0,0,1],[0,1,1],...
                            [1,0,1],[0.8500 0.3250 0.0980],[0,1,0]};
            var_list={'S_0','S_hat','w_hat','u_tilde'};
            xlabel_list={'$\bar{S}_0+z$','$\widehat{S}$','$\widehat{w}$','$\widetilde{u}$'};
            for var_ind=1:length(var_list)
                var_name=var_list{var_ind};
                for slurm_ind=1:length(slurm_num)
                    if strcmp(var_list{var_ind},'S_0')
                        data{slurm_ind}.x=dedalus_post_my{slurm_ind,1}.(var_name)+dedalus_post_my{slurm_ind,1}.z_list; 
                    else
                        data{slurm_ind}.x=dedalus_post_my{slurm_ind,1}.(var_name); 
                    end
                    data{slurm_ind}.y=dedalus_post_my{slurm_ind,1}.z_list;
                end
                data=data(order);
                data_sym=data(1:4);
                data_asym=data(5:end);
                plot_config.label_list={1,xlabel_list{var_ind},'$z$'};
                plot_config.print_size=[1,600,1100];
                plot_config.name=['C:\Figure\DDC_LST\',group_name,'_profile_sym_',var_name,'.png'];
                plot_config.RGB={1,RGB{1:4}};
                plot_config.legend_list={1,legend_list{1:4}};
                plot_config.fontsize_legend=37;
                plot_config.linewidth=5;
                if var_ind>1 %%get rid of y label and the legend
                    plot_config.legend_list={0}; plot_config.label_list{3}='';
                end
                plot_line(data_sym,plot_config);
                
                plot_config.RGB={1,RGB{5:end}};
                plot_config.legend_list={1,legend_list{5:end}};
                if var_ind>1 %get rid of the y label and legend for the following results.
                    plot_config.legend_list={0}; plot_config.label_list{3}='';
                end
                plot_config.name=['C:\Figure\DDC_LST\',group_name,'_profile_asym_',var_name,'.png'];
                plot_line(data_asym,plot_config);
            end
            
        case 'HB_benard_salt_finger_profile_S1_yang'
            z_list=dedalus_post_my{1,1}.z_list;
            slurm_len=length(slurm_num);
            for slurm_ind=1:length(slurm_num)
                for content_ind=1:size(dedalus_post_my,2)
                    if ~isempty(dedalus_post_my{slurm_ind,content_ind})
                        R_rho_T2S(content_ind+(slurm_ind-1)*slurm_len)=dedalus_post_my{slurm_ind,content_ind}.Ra_T/dedalus_post_my{slurm_ind,content_ind}.Ra_S2T;
                        S_0(:,content_ind+(slurm_ind-1)*slurm_len)=dedalus_post_my{slurm_ind,content_ind}.S_0+z_list;
                        T_0(:,content_ind+(slurm_ind-1)*slurm_len)=dedalus_post_my{slurm_ind,content_ind}.T_0+z_list;
                    end
                end
            end
%             R_rho_T2S_list=[0.1,0.2,0.5,1,2,5,10];
            R_rho_T2S_list=[10,5,2,1,0.5,0.2,0.1];
            for R_rho_T2S_ind=1:length(R_rho_T2S_list)
               ind=find(abs(R_rho_T2S-R_rho_T2S_list(R_rho_T2S_ind))<0.01)
               data{R_rho_T2S_ind}.x=S_0(:,ind(1))+R_rho_T2S_ind-1;
               data{R_rho_T2S_ind}.y=z_list;
               data{R_rho_T2S_ind+length(R_rho_T2S_list)}.x=T_0(:,ind(1))+R_rho_T2S_ind-1;
               data{R_rho_T2S_ind+length(R_rho_T2S_list)}.y=z_list;
            end

            plot_config.label_list={1,'','$z$'};
            plot_config.Markerindex=3;
            plot_config.user_color_style_marker_list={'k-','k-','k-','k-','k-','k-','k-','b--','b--','b--','b--','b--','b--','b--'};
            plot_config.linewidth=3;
            plot_config.name=['C:\Figure\DDC_LST\',group_name,'_Ra_1e5.png'];
            plot_config.xlim_list=[1,-0.1,7.1];
            plot_config.ylim_list=[1,-0.05,1.05];
            plot_config.xtick_list=[1,0,1,2,3,4,5,6,7];
            plot_config.ytick_list=[1,0,0.5,1];
            plot_config.print_size=[1,1500,800];
            plot_line(data,plot_config);

        case 'HB_benard_salt_finger_R_rho_T2S'
            for slurm_ind=1:size(dedalus_post_my,1)
                for content_ind=1:size(dedalus_post_my,2)
                    if ~isempty(dedalus_post_my{slurm_ind,content_ind})
                        data{slurm_ind}.x(content_ind)=dedalus_post_my{slurm_ind,content_ind}.Ra_T/dedalus_post_my{slurm_ind,content_ind}.Ra_S2T;
                        %R_rho_T2S_list_content(content_ind)=data{slurm_ind}.x(content_ind);
                        data{slurm_ind}.y(content_ind)=abs(dedalus_post_my{slurm_ind,content_ind}.Nu_S(1));
                        
                        eig_inf_ind=find(real(dedalus_post_my{slurm_ind,content_ind}.eigenvalues)>10^6);
                        dedalus_post_my{slurm_ind,content_ind}.eigenvalues(eig_inf_ind)=[];
                        lambda_max{slurm_ind}.y(content_ind)=dedalus_post_my{slurm_ind,content_ind}.eigenvalues(1);
                   
                    end
                end
                %[~,ind]=sort(data{slurm_ind}.x);
                %dedalus_post_my(slurm_ind,:)=dedalus_post_my(slurm_ind,ind);
            end
            sol_num=7;
            for sol_ind=1:sol_num
                %%combine the data that increase R_rho_S2T and decrease
                %%R_rho_S2T
                data_all{sol_ind}.x=[data{sol_ind}.x,data{sol_ind+sol_num}.x,data{sol_ind+2*sol_num}.x];
                R_rho_T2S_list_content=data_all{sol_ind}.x;
                data_all{sol_ind}.y=[data{sol_ind}.y,data{sol_ind+sol_num}.y,data{sol_ind+2*sol_num}.y];
                lambda_max_all{sol_ind}.y=[lambda_max{sol_ind}.y,lambda_max{sol_ind+sol_num}.y,lambda_max{sol_ind+2*sol_num}.y];
                
                %sort data
                [~,ind]=sort(data_all{sol_ind}.x);
%                 ind=ind(2:end);
                data_all{sol_ind}.x=data_all{sol_ind}.x(ind);
                data_all{sol_ind}.y=data_all{sol_ind}.y(ind);
                lambda_max_all{sol_ind}.x=data_all{sol_ind}.x;
                lambda_max_all{sol_ind}.y=lambda_max_all{sol_ind}.y(ind);
                stable_ind=find(real(lambda_max_all{sol_ind}.y)<=0);
                if isempty(stable_ind)
                    data_all_stable{sol_ind}.x=NaN;
                    data_all_stable{sol_ind}.y=NaN;
                    data_all_stable_2{sol_ind}.x=NaN;
                    data_all_stable_2{sol_ind}.y=NaN;
                elseif isempty(find(diff(stable_ind)>1))
                    data_all_stable{sol_ind}.x=data_all{sol_ind}.x(stable_ind);
                    data_all_stable{sol_ind}.y=data_all{sol_ind}.y(stable_ind);
                    data_all_stable_2{sol_ind}.x=NaN; data_all_stable_2{sol_ind}.y=NaN;
                else
                    stable_2_ind=find(diff(stable_ind)>1);
                    stable_2_ind=stable_2_ind(1);
                    data_all_stable{sol_ind}.x=data_all{sol_ind}.x(stable_ind(1:stable_2_ind));
                    data_all_stable{sol_ind}.y=data_all{sol_ind}.y(stable_ind(1:stable_2_ind));
                    data_all_stable_2{sol_ind}.x=data_all{sol_ind}.x(stable_ind(1+stable_2_ind:end));
                    data_all_stable_2{sol_ind}.y=data_all{sol_ind}.y(stable_ind(1+stable_2_ind:end));
                end
                
                unstable_ind=find(real(lambda_max_all{sol_ind}.y)>0);
                if isempty(unstable_ind)
                    data_all_unstable{sol_ind}.x=NaN;
                    data_all_unstable{sol_ind}.y=NaN;
                    data_all_unstable_2{sol_ind}.x=NaN; data_all_unstable_2{sol_ind}.y=NaN;
                elseif isempty(find(diff(unstable_ind)>1))
                    data_all_unstable{sol_ind}.x=data_all{sol_ind}.x(unstable_ind);
                    data_all_unstable{sol_ind}.y=data_all{sol_ind}.y(unstable_ind);
                    data_all_unstable_2{sol_ind}.x=NaN; data_all_unstable_2{sol_ind}.y=NaN;
                else
                    unstable_2_ind=find(diff(unstable_ind)>1);
                    unstable_2_ind=unstable_2_ind(1);
                    data_all_unstable{sol_ind}.x=data_all{sol_ind}.x(unstable_ind(1:unstable_2_ind));
                    data_all_unstable{sol_ind}.y=data_all{sol_ind}.y(unstable_ind(1:unstable_2_ind));
                    data_all_unstable_2{sol_ind}.x=data_all{sol_ind}.x(unstable_ind(1+unstable_2_ind:end));
                    data_all_unstable_2{sol_ind}.y=data_all{sol_ind}.y(unstable_ind(1+unstable_2_ind:end));
                end
%                 if all(real(lambda_max{sol_ind}.y)>0)
%                     %plot_config.user_color_style_marker_list{sol_ind}=[color_list{sol_ind},'--'];
%                 elseif all(real(lambda_max{sol_ind}.y)<0)
%                     %plot_config.user_color_style_marker_list{sol_ind}=[color_list{sol_ind},'-'];
%                 else 
%                     error('There is a transition of stability in parameter regime');
%                 end
            end
            Nu_salt_finger_yang=dedalus_post_my{1,1}.get_Nu_salt_finger_yang();
            plot_config.legend_list={1,'S2','A2','S2L','A3','S1','A1','S3'};
            order=[5,1,7,3,6,2,4];%This will give S1, S2, S3, S2L, A1, A2, A3 order...
            plot_config.legend_list(2:end)=plot_config.legend_list(order+1);
            data_R_rho_T2S=[data_all_stable(order),data_all_stable_2(order),data_all_unstable(order),data_all_unstable_2(order)];
            plot_config.user_color_style_marker_list={'-','-','-','-','-','-','-',...
                                                    '-','-','-','-','-','-','-',...                                        
                                                    ':',':',':',':',':',':',':',...
                                                    ':',':',':',':',':',':',':'};
            data_len=length(data_R_rho_T2S);
            data_R_rho_T2S{data_len+1}.x=Nu_salt_finger_yang.Ra_T_1e5_Nu_S(:,1);
            data_R_rho_T2S{data_len+1}.y=Nu_salt_finger_yang.Ra_T_1e5_Nu_S(:,2);
            plot_config.user_color_style_marker_list{data_len+1}='k*';
            plot_config.Markerindex=3;
            plot_config.label_list={1,'$R_\rho$','$Sh$'};
            plot_config.loglog=[1,1];
            plot_config.ylim_list=[1,1,100];
            plot_config.ytick_list=[1,1,2,5,10,20,50,100];
            plot_config.xlim_list=[1,0.2,500];
            plot_config.xtick_list=[1,0.1,1,10,100,500];
            plot_config.linewidth=5;
            plot_config.RGB={1,[0,0,0],[1,0,0],[0,0,1],[0,1,1],[1,0,1],...
                                [0.8500 0.3250 0.0980],[0,1,0],...
                                [0,0,0],[1,0,0],[0,0,1],[0,1,1],[1,0,1],...
                                [0.8500 0.3250 0.0980],[0,1,0],...
                                [0,0,0],[1,0,0],[0,0,1],[0,1,1],[1,0,1],...
                                [0.8500 0.3250 0.0980],[0,1,0],...
                                [0,0,0],[1,0,0],[0,0,1],[0,1,1],[1,0,1],...
                                [0.8500 0.3250 0.0980],[0,1,0],...
                                [0,0,0]};
            plot_config.name=['C:\Figure\DDC_LST\',group_name,'_Sh_R_rho.png'];
            plot_line(data_R_rho_T2S,plot_config);
            
            %plot the profile where the staircase solution gain stability
            order=[5,1,7,3,6,2,4];%This will give S1, S2, S3, S2L, A1, A2, A3 order...
%             R_rho_T2S_list=[30,40];
%             for content_ind=1:size(dedalus_post_my,2)
%                 R_rho_T2S_list_content=dedalus_post_my{1,content_ind}.;
%             end
            content_list=[16,17];%for R_rho_T2S=30 and 40
            clear data;
            plot_config.ylim_list=[1,0,1];
            plot_config.xlim_list=0;
            plot_config.xtick_list=0; plot_config.ytick_list=0
            plot_config.loglog=[0,0];
            plot_config.Markerindex=1;
            for content_ind=content_list
                %R_rho_T2S=R_rho_T2S_list(R_rho_T2S_ind);
                %content_ind=find(abs(R_rho_T2S-R_rho_T2S_list_content)<0.01);
                var_list={'S_0','S_hat','w_hat','u_tilde'};
                xlabel_list={'$\bar{S}_0+z$','$\widehat{S}$','$\widehat{w}$','$\widetilde{u}$'};
                R_rho_T2S=dedalus_post_my{14+1,content_ind}.Ra_T/dedalus_post_my{14+1,content_ind}.Ra_S2T;
                plot_config.legend_list={1,'S1','S2','S3','S2L'};
                for var_ind=1:length(var_list)
                    var_name=var_list{var_ind};
                    for slurm_ind=1:4
                        if strcmp(var_list{var_ind},'S_0')
                            data{slurm_ind}.x=dedalus_post_my{14+order(slurm_ind),content_ind}.(var_name)+dedalus_post_my{14+order(slurm_ind),1}.z_list; 
                        else
                            data{slurm_ind}.x=dedalus_post_my{14+order(slurm_ind),content_ind}.(var_name); 
                        end
                        data{slurm_ind}.y=dedalus_post_my{slurm_ind,1}.z_list;
                    end
%                     data=data(order);
                    data_sym=data(1:4);
%                     data_asym=data(5:end);
                    plot_config.label_list={1,xlabel_list{var_ind},'$z$'};
                    plot_config.print_size=[1,600,1100];
                    plot_config.name=['C:\Figure\DDC_LST\',group_name,'_profile_sym_R_rho_S2T',num2str(R_rho_T2S),'_',var_name,'.png'];
                    plot_config.fontsize_legend=37;
                    plot_config.linewidth=5;
                    if var_ind>1 %%get rid of y label and the legend
                        plot_config.legend_list={0}; plot_config.label_list{3}='';
                    end
                    plot_line(data_sym,plot_config);

                end
            end
            
            
        case 'HB_benard_salt_finger_tau'
            for slurm_ind=1:size(dedalus_post_my,1)
                for content_ind=1:size(dedalus_post_my,2)
                    if ~isempty(dedalus_post_my{slurm_ind,content_ind})
                        data{slurm_ind}.x(content_ind)=dedalus_post_my{slurm_ind,content_ind}.tau;
                        data{slurm_ind}.y(content_ind)=abs(dedalus_post_my{slurm_ind,content_ind}.Nu_S(1));
                        
                        %Update 2022/03/02, remove all eigenvalues that are
                        %particularly large. This is because we solve an
                        %generalized eigenvalue problem, and some inf
                        %eigenvalue may be not filtterd.... I should also
                        %modify the eigenvectors...
                        eig_inf_ind=find(real(dedalus_post_my{slurm_ind,content_ind}.eigenvalues)>10^6);
                        dedalus_post_my{slurm_ind,content_ind}.eigenvalues(eig_inf_ind)=[];
                        lambda_max{slurm_ind}.y(content_ind)=dedalus_post_my{slurm_ind,content_ind}.eigenvalues(1);
                    end
                end
                %[~,ind]=sort(data{slurm_ind,:}.x);
                %dedalus_post_my(slurm_ind,:)=dedalus_post_my(slurm_ind,ind);
            end
            sol_num=7;
            for sol_ind=1:sol_num
                %%combine the data that increase R_rho_S2T and decrease
                %%R_rho_S2T
                %data_all{sol_ind}.x=[data{sol_ind}.x,data{sol_ind+sol_num}.x,data{sol_ind+2*sol_num}.x];
                %data_all{sol_ind}.y=[data{sol_ind}.y,data{sol_ind+sol_num}.y,data{sol_ind+2*sol_num}.y];
                data_all{sol_ind}=data{sol_ind};
                lambda_max_all{sol_ind}=lambda_max{sol_ind};
                %lambda_max{sol_ind}.y=[lambda_max{sol_ind}.y,lambda_max{sol_ind+sol_num}.y,lambda_max{sol_ind+sol_num}.y,data{sol_ind+2*sol_num}.y];
                
                %sort data
                [~,ind]=sort(data_all{sol_ind}.x);
%                 ind=ind(2:end);
                data_all{sol_ind}.x=data_all{sol_ind}.x(ind);
                data_all{sol_ind}.y=data_all{sol_ind}.y(ind);
                lambda_max_all{sol_ind}.x=data_all{sol_ind}.x;
                lambda_max_all{sol_ind}.y=lambda_max_all{sol_ind}.y(ind);
                unstable_ind=find(real(lambda_max_all{sol_ind}.y)>0);
                if isempty(unstable_ind)
                    data_all_unstable{sol_ind}.x=NaN;
                    data_all_unstable{sol_ind}.y=NaN;
                    data_all_unstable_2{sol_ind}.x=NaN;
                    data_all_unstable_2{sol_ind}.y=NaN;
                elseif isempty(find(diff(unstable_ind)>1))
                    data_all_unstable{sol_ind}.x=data_all{sol_ind}.x(unstable_ind);
                    data_all_unstable{sol_ind}.y=data_all{sol_ind}.y(unstable_ind);
                    data_all_unstable_2{sol_ind}.x=NaN; data_all_unstable_2{sol_ind}.y=NaN;
                else
                    unstable_2_ind=find(diff(unstable_ind)>1);
                    unstable_2_ind=unstable_2_ind(1);
                    data_all_unstable{sol_ind}.x=data_all{sol_ind}.x(unstable_ind(1:unstable_2_ind));
                    data_all_unstable{sol_ind}.y=data_all{sol_ind}.y(unstable_ind(1:unstable_2_ind));
                    data_all_unstable_2{sol_ind}.x=data_all{sol_ind}.x(unstable_ind(1+unstable_2_ind:end));
                    data_all_unstable_2{sol_ind}.y=data_all{sol_ind}.y(unstable_ind(1+unstable_2_ind:end));
                end
                
                stable_ind=find(real(lambda_max_all{sol_ind}.y<=0));
                if isempty(stable_ind)
                    data_all_stable{sol_ind}.x=NaN;
                    data_all_stable{sol_ind}.y=NaN;
                    data_all_stable_2{sol_ind}.x=NaN;
                    data_all_stable_2{sol_ind}.y=NaN;
                elseif isempty(find(diff(stable_ind)>1))
                    data_all_stable{sol_ind}.x=data_all{sol_ind}.x(stable_ind);
                    data_all_stable{sol_ind}.y=data_all{sol_ind}.y(stable_ind);
                    data_all_stable_2{sol_ind}.x=NaN; data_all_stable_2{sol_ind}.y=NaN;
                else
                    stable_2_ind=find(diff(stable_ind)>1);
                    stable_2_ind=stable_2_ind(1);
                    data_all_stable{sol_ind}.x=data_all{sol_ind}.x(stable_ind(1:stable_2_ind));
                    data_all_stable{sol_ind}.y=data_all{sol_ind}.y(stable_ind(1:stable_2_ind));
                    data_all_stable_2{sol_ind}.x=data_all{sol_ind}.x(stable_ind(1+stable_2_ind:end));
                    data_all_stable_2{sol_ind}.y=data_all{sol_ind}.y(stable_ind(1+stable_2_ind:end));
                end
%                 if all(real(lambda_max{sol_ind}.y)>0)
%                     %plot_config.user_color_style_marker_list{sol_ind}=[color_list{sol_ind},'--'];
%                 elseif all(real(lambda_max{sol_ind}.y)<0)
%                     %plot_config.user_color_style_marker_list{sol_ind}=[color_list{sol_ind},'-'];
%                 else 
%                     error('There is a transition of stability in parameter regime');
%                 end
            end
            %Nu_salt_finger_yang=dedalus_post_my{1,1}.get_Nu_salt_finger_yang();
            plot_config.legend_list={1,'S2','A2','S2L','A3','S1','A1','S3'};
            order=[5,1,7,3,6,2,4];%This will give S1, S2, S3, S2L, A1, A2, A3 order...
            plot_config.legend_list(2:end)=plot_config.legend_list(order+1);
            data_tau=[data_all_stable(order),data_all_stable_2(order),data_all_unstable(order),data_all_unstable_2(order)];
            plot_config.user_color_style_marker_list={'-','-','-','-','-','-','-',...
                                                        '-','-','-','-','-','-','-',...                                          
                                                        ':',':',':',':',':',':',':',...
                                                      ':',':',':',':',':',':',':'};
            %data_len=length(data_R_rho_T2S);
            %data_R_rho_T2S{data_len+1}.x=Nu_salt_finger_yang.Ra_T_1e5_Nu_S(:,1);
            %data_R_rho_T2S{data_len+1}.y=Nu_salt_finger_yang.Ra_T_1e5_Nu_S(:,2);
            %plot_config.user_color_style_marker_list{data_len+1}='k*';
            plot_config.Markerindex=3;
            plot_config.label_list={1,'$\tau$','$Sh$'};
            plot_config.loglog=[1,1];
            plot_config.ylim_list=[1,1,100];
            plot_config.ytick_list=[1,1,2,5,10,20,50,100];
            plot_config.xlim_list=[1,0.01,1];
            plot_config.xtick_list=[1,0.01,0.1,0.3,0.5,0.8,1];
            plot_config.linewidth=5;
            plot_config.RGB={1,[0,0,0],[1,0,0],[0,0,1],[0,1,1],[1,0,1],...
                                [0.8500 0.3250 0.0980],[0,1,0],...
                                [0,0,0],[1,0,0],[0,0,1],[0,1,1],[1,0,1],...
                                [0.8500 0.3250 0.0980],[0,1,0],...
                                [0,0,0],[1,0,0],[0,0,1],[0,1,1],[1,0,1],...
                                [0.8500 0.3250 0.0980],[0,1,0],...
                                [0,0,0],[1,0,0],[0,0,1],[0,1,1],[1,0,1],...
                                [0.8500 0.3250 0.0980],[0,1,0],...
                                [0,0,0]};
            plot_config.name=['C:\Figure\DDC_LST\',group_name,'_Sh_tau.png'];
            plot_line(data_tau,plot_config);
            
            
            %----------
            %plot the profile where the staircase solution gain stability
%             R_rho_T2S_list=[30,40];
%             for content_ind=1:size(dedalus_post_my,2)
%                 R_rho_T2S_list_content=dedalus_post_my{1,content_ind}.;
%             end
%             content_list=[11];%for R_rho_T2S=30 and 40
            
            clear data;
%             plot_config.legend_list=plot_config.legend_list(1:5);
            plot_config.ylim_list=[1,0,1];
            plot_config.xlim_list=0;
            plot_config.xtick_list=0; plot_config.ytick_list=0
            plot_config.loglog=[0,0];
            plot_config.Markerindex=1;
            
            tau=0.18;
            order=[5,1,7,3,6,2,4];%This will give S1, S2, S3, S2L, A1, A2, A3 order...
            content_list=[11,11,10,11];%This make sure dedalus_post_my{order(slurm_ind),content_list(slurm_ind)} gives the right tau=0.18
            %R_rho_T2S=R_rho_T2S_list(R_rho_T2S_ind);
            %content_ind=find(abs(R_rho_T2S-R_rho_T2S_list_content)<0.01);
            var_list={'S_0','S_hat','w_hat','u_tilde'};
            xlabel_list={'$\bar{S}_0+z$','$\widehat{S}$','$\widehat{w}$','$\widetilde{u}$'};
            %R_rho_T2S=dedalus_post_my{14+1,content_ind}.Ra_T/dedalus_post_my{14+1,content_ind}.Ra_S2T;
            plot_config.legend_list={1,'S1','S2','S3','S2L'};
            for var_ind=1:length(var_list)
                var_name=var_list{var_ind};
                for slurm_ind=1:4
                    if strcmp(var_list{var_ind},'S_0')
                        data{slurm_ind}.x=dedalus_post_my{order(slurm_ind),content_list(slurm_ind)}.(var_name)+dedalus_post_my{order(slurm_ind),1}.z_list; 
                    else
                        data{slurm_ind}.x=dedalus_post_my{order(slurm_ind),content_list(slurm_ind)}.(var_name); 
                    end
                    data{slurm_ind}.y=dedalus_post_my{slurm_ind,1}.z_list;
                end
%                     data=data(order);
                data_sym=data(1:4);
%                     data_asym=data(5:end);
                plot_config.label_list={1,xlabel_list{var_ind},'$z$'};
                plot_config.print_size=[1,600,1100];
                plot_config.name=['C:\Figure\DDC_LST\',group_name,'_profile_sym_tau',num2str(tau),'_',var_name,'.png'];
                plot_config.fontsize_legend=37;
                plot_config.linewidth=5;
                if var_ind>1 %%get rid of y label and the legend
                    plot_config.legend_list={0}; plot_config.label_list{3}='';
                end
                plot_line(data_sym,plot_config);

            end

            
                
            
        case 'trevisan_tau'
            for slurm_ind=1:size(dedalus_post_my,1)
                for content_ind=1:size(dedalus_post_my,2)
                    tau(slurm_ind,content_ind)=dedalus_post_my{slurm_ind,content_ind}.tau;
                    Nu_S(slurm_ind,content_ind)=dedalus_post_my{slurm_ind,content_ind}.Nu_S(1);
                    
                end
                [tau_tmp,ind]=sort(tau(slurm_ind,:));
                tau(slurm_ind,:)=tau_tmp;  
                Nu_S(slurm_ind,:)=Nu_S(slurm_ind,ind);
                data{slurm_ind}.x=1./tau(slurm_ind,:);
                data{slurm_ind}.y=Nu_S(slurm_ind,:);
                eta_Nu_S_Le(slurm_ind)=scaling(data{slurm_ind}.x(1:4),data{slurm_ind}.y(1:4));
            end
            porous_trevisan_tau=dedalus_post_my{1,1}.get_porous_trevisan_tau
            data{6}.x=porous_trevisan_tau.Sh_Le_Ra_50(:,1);
            data{6}.y=porous_trevisan_tau.Sh_Le_Ra_50(:,2);
            data{7}.x=porous_trevisan_tau.Sh_Le_Ra_100(:,1);
            data{7}.y=porous_trevisan_tau.Sh_Le_Ra_100(:,2);
            data{8}.x=porous_trevisan_tau.Sh_Le_Ra_200(:,1);
            data{8}.y=porous_trevisan_tau.Sh_Le_Ra_200(:,2);
            data{9}.x=porous_trevisan_tau.Sh_Le_Ra_400(:,1);
            data{9}.y=porous_trevisan_tau.Sh_Le_Ra_400(:,2);
            data{10}.x=porous_trevisan_tau.Sh_Le_Ra_1000(:,1);
            data{10}.y=porous_trevisan_tau.Sh_Le_Ra_1000(:,2);
            plot_config.loglog=[1,1];
            plot_config.Markerindex=3;
            plot_config.user_color_style_marker_list={'k-','b--','r:','m-','g--','ko','bsquare','r^','m*','gdiamond'};
            plot_config.label_list={1,'Le','Sh'};
            plot_config.name=['C:\Figure\DDC_LST\',group_name,'Sh_Le.png'];
            plot_config.print_size=[1,1000,900];
            plot_config.legend_list={1,'Ra=50','Ra=100','Ra=200','Ra=400','Ra=1000','Ra=50 (DNS)','Ra=100 (DNS)','Ra=200 (DNS)','Ra=400 (DNS)','Ra=1000 (DNS)'}
            plot_config.fontsize_legend=22;
            plot_line(data,plot_config);
        case 'trevisan_contour'
            
        case 'rosenberg_tau'
            for slurm_ind=1:size(dedalus_post_my,1)
                for content_ind=1:size(dedalus_post_my,2)
                    tau(slurm_ind,content_ind)=dedalus_post_my{slurm_ind,content_ind}.tau;
                    Nu(slurm_ind,content_ind)=dedalus_post_my{slurm_ind,content_ind}.Nu(1);
                    Nu_S(slurm_ind,content_ind)=dedalus_post_my{slurm_ind,content_ind}.Nu_S(1);
                end
                [tau_tmp,ind]=sort(tau(slurm_ind,:));
                tau(slurm_ind,:)=tau_tmp;  
                Nu_S(slurm_ind,:)=Nu_S(slurm_ind,ind);
                Nu(slurm_ind,:)=Nu(slurm_ind,ind);
                data_Nu_S{slurm_ind}.x=1./tau(slurm_ind,:);
                data_Nu_S{slurm_ind}.y=Nu_S(slurm_ind,:);
                data_Nu{slurm_ind}.x=1./tau(slurm_ind,:);
                data_Nu{slurm_ind}.y=Nu(slurm_ind,:);
                eta_Nu_S_Le(slurm_ind)=scaling(data_Nu_S{slurm_ind}.x(1:4),data_Nu_S{slurm_ind}.y(1:4));
            end
            porous_rosenberg_tau=dedalus_post_my{1,1}.get_porous_rosenberg_tau()
            data_Nu{5}.x=porous_rosenberg_tau.Nu_Le_Ra_100(:,1);
            data_Nu{5}.y=porous_rosenberg_tau.Nu_Le_Ra_100(:,2);
            data_Nu{6}.x=porous_rosenberg_tau.Nu_Le_Ra_150(:,1);
            data_Nu{6}.y=porous_rosenberg_tau.Nu_Le_Ra_150(:,2);
            data_Nu{7}.x=porous_rosenberg_tau.Nu_Le_Ra_300(:,1);
            data_Nu{7}.y=porous_rosenberg_tau.Nu_Le_Ra_300(:,2);
            data_Nu{8}.x=porous_rosenberg_tau.Nu_Le_Ra_600(:,1);
            data_Nu{8}.y=porous_rosenberg_tau.Nu_Le_Ra_600(:,2);
         
            data_Nu_S{5}.x=porous_rosenberg_tau.Sh_Le_Ra_100(:,1);
            data_Nu_S{5}.y=porous_rosenberg_tau.Sh_Le_Ra_100(:,2);
            data_Nu_S{6}.x=porous_rosenberg_tau.Sh_Le_Ra_150(:,1);
            data_Nu_S{6}.y=porous_rosenberg_tau.Sh_Le_Ra_150(:,2);
            data_Nu_S{7}.x=porous_rosenberg_tau.Sh_Le_Ra_300(:,1);
            data_Nu_S{7}.y=porous_rosenberg_tau.Sh_Le_Ra_300(:,2);
            data_Nu_S{8}.x=porous_rosenberg_tau.Sh_Le_Ra_600(:,1);
            data_Nu_S{8}.y=porous_rosenberg_tau.Sh_Le_Ra_600(:,2);
            plot_config.loglog=[1,1];
            plot_config.Markerindex=3;
            plot_config.user_color_style_marker_list={'k-','b--','r:','m-','ko','bsquare','r^','m*'};
            plot_config.label_list={1,'Le','Nu'};
            plot_config.print_size=[1,1000,900];
            plot_config.legend_list={1,'Ra=100','Ra=150','Ra=300','Ra=600','Ra=100 (DNS)','Ra=150 (DNS)','Ra=300 (DNS)','Ra=600 (DNS)'};
            plot_config.fontsize_legend=22;
            plot_config.loglog=[1,1];
            plot_config.xlim_list=[1,10,100];
            plot_config.ylim_list=[1,1,8];
            plot_config.name=['C:\Figure\DDC_LST\',group_name,'Nu_Le.png'];
            plot_config.fontsize_legend=14;
            plot_line(data_Nu,plot_config);
            
            plot_config.xlim_list=[1,9.9,101];
            plot_config.xtick_list=[1,10,100];
            plot_config.ylim_list=[1,10,100];
            plot_config.loglog=[1,1];
            plot_config.label_list={1,'Le','Sh'};
            plot_config.name=['C:\Figure\DDC_LST\',group_name,'Sh_Le.png'];
            plot_config.fontsize_legend=18;
            plot_line(data_Nu_S,plot_config);
        
        case 'rosenberg_R_rho'
            for slurm_ind=1:size(dedalus_post_my,1)
                for content_ind=1:size(dedalus_post_my,2)
                    R_rho_S2T(slurm_ind,content_ind)=dedalus_post_my{slurm_ind,content_ind}.Ra_S2T/dedalus_post_my{slurm_ind,content_ind}.Ra_T;
                    Nu(slurm_ind,content_ind)=dedalus_post_my{slurm_ind,content_ind}.Nu(1);
                    Nu_S(slurm_ind,content_ind)=dedalus_post_my{slurm_ind,content_ind}.Nu_S(1);
                end
                [R_rho_S2T_tmp,ind]=sort(R_rho_S2T(slurm_ind,:));
                R_rho_S2T(slurm_ind,:)=R_rho_S2T_tmp;  
                Nu_S(slurm_ind,:)=Nu_S(slurm_ind,ind);
                Nu(slurm_ind,:)=Nu(slurm_ind,ind);
                data_Nu_S{slurm_ind}.x=R_rho_S2T(slurm_ind,:);
                data_Nu_S{slurm_ind}.y=Nu_S(slurm_ind,:);
                data_Nu{slurm_ind}.x=R_rho_S2T(slurm_ind,:);
                data_Nu{slurm_ind}.y=Nu(slurm_ind,:);
            end
            porous_rosenberg_R_rho_S2T=dedalus_post_my{1,1}.get_porous_rosenberg_R_rho_S2T();
            data_Nu{5}.x=porous_rosenberg_R_rho_S2T.Nu_R_rho_Ra_100(:,1);
            data_Nu{5}.y=porous_rosenberg_R_rho_S2T.Nu_R_rho_Ra_100(:,2);
            data_Nu{6}.x=porous_rosenberg_R_rho_S2T.Nu_R_rho_Ra_150(:,1);
            data_Nu{6}.y=porous_rosenberg_R_rho_S2T.Nu_R_rho_Ra_150(:,2);
            data_Nu{7}.x=porous_rosenberg_R_rho_S2T.Nu_R_rho_Ra_300(:,1);
            data_Nu{7}.y=porous_rosenberg_R_rho_S2T.Nu_R_rho_Ra_300(:,2);
            data_Nu{8}.x=porous_rosenberg_R_rho_S2T.Nu_R_rho_Ra_600(:,1);
            data_Nu{8}.y=porous_rosenberg_R_rho_S2T.Nu_R_rho_Ra_600(:,2);
         
            data_Nu_S{5}.x=porous_rosenberg_R_rho_S2T.Sh_R_rho_Ra_100(:,1);
            data_Nu_S{5}.y=porous_rosenberg_R_rho_S2T.Sh_R_rho_Ra_100(:,2);
            data_Nu_S{6}.x=porous_rosenberg_R_rho_S2T.Sh_R_rho_Ra_150(:,1);
            data_Nu_S{6}.y=porous_rosenberg_R_rho_S2T.Sh_R_rho_Ra_150(:,2);
            data_Nu_S{7}.x=porous_rosenberg_R_rho_S2T.Sh_R_rho_Ra_300(:,1);
            data_Nu_S{7}.y=porous_rosenberg_R_rho_S2T.Sh_R_rho_Ra_300(:,2);
            data_Nu_S{8}.x=porous_rosenberg_R_rho_S2T.Sh_R_rho_Ra_600(:,1);
            data_Nu_S{8}.y=porous_rosenberg_R_rho_S2T.Sh_R_rho_Ra_600(:,2);
            plot_config.loglog=[1,1];
            plot_config.Markerindex=3;
            plot_config.user_color_style_marker_list={'k-','b--','r:','m-','ko','bsquare','r^','m*'};
            plot_config.label_list={1,'$R_\rho^*$','Nu'};
            plot_config.print_size=[1,1000,900];
            plot_config.legend_list={1,'Ra=100','Ra=150','Ra=300','Ra=600','Ra=100 (DNS)','Ra=150 (DNS)','Ra=300 (DNS)','Ra=600 (DNS)'};
            plot_config.fontsize_legend=22;
            plot_config.loglog=[1,1];
            plot_config.xlim_list=[1,0,0.401];
            plot_config.ylim_list=[1,0,15];
            plot_config.name=['C:\Figure\DDC_LST\',group_name,'Nu_R_rho.png'];
            plot_config.fontsize_legend=20;
            plot_config.loglog=[0,0];
            plot_line(data_Nu,plot_config);
            
            plot_config.xlim_list=[1,0,0.401];
            plot_config.ylim_list=[1,0,80];
            plot_config.loglog=[1,1];
            plot_config.label_list={1,'$R_\rho^*$','Sh'};
            plot_config.name=['C:\Figure\DDC_LST\',group_name,'Sh_R_rho.png'];
            plot_config.fontsize_legend=20;
            plot_config.loglog=[0,0];
            plot_line(data_Nu_S,plot_config);
            
            
        case 'test'
            for content_ind=1:size(dedalus_post_my,2)
                data_Nu{1}.x(content_ind)=dedalus_post_my{content_ind}.HB_porous_3_layer_h;
                data_Nu{1}.y(content_ind)=dedalus_post_my{content_ind}.Nu(1);
            end
    end

end
% 
% for slurm_ind=length(slurm_num)
%     h5_name=['C:\Data\dedalus\dedalus_',...
%             slurm_num{slurm_ind},...
%                 '\analysis\analysis_s1.h5'];
% %             '\'];
% %     files = dir(h5_name);
% %     for folder_ind=1:length(files)
% %         h5_name=
%             
%     %     h5_name=['C:\Data\dedalus\dedalus_',...
%     %             slurm_num{slurm_ind},...
%     %             '\data.h5'];
%         h5disp(h5_name);
%         obj.z_list=h5read(h5_name,'/scales/z/1.0');
%         obj.w_hat=h5read(h5_name,'/tasks/w_hat');
%         obj.p_hat=h5read(h5_name,'/tasks/p_hat');
%         obj.T_hat=h5read(h5_name,'/tasks/T_hat');
%         obj.d_T_hat=h5read(h5_name,'/tasks/d_T_hat');
%         obj.S_hat=h5read(h5_name,'/tasks/S_hat');
%         obj.d_S_hat=h5read(h5_name,'/tasks/d_S_hat');
%         try
%             obj.w_hat_2=h5read(h5_name,'/tasks/w_hat_2');
%             obj.p_hat_2=h5read(h5_name,'/tasks/p_hat_2');
%             obj.T_hat_2=h5read(h5_name,'/tasks/T_hat_2');
%             obj.d_T_hat_2=h5read(h5_name,'/tasks/d_T_hat_2');
%             obj.S_hat_2=h5read(h5_name,'/tasks/S_hat_2');
%             obj.d_S_hat_2=h5read(h5_name,'/tasks/d_S_hat_2');
%         catch
%             disp('No second harmonic')
%         end
%         obj.T_0=h5read(h5_name,'/tasks/T_0');
%         obj.d_T_0=h5read(h5_name,'/tasks/d_T_0');
%         obj.S_0=h5read(h5_name,'/tasks/S_0');
%         obj.d_S_0=h5read(h5_name,'/tasks/d_S_0');
% 
%         obj.dy_T_mean=-1;
%         obj.dy_S_mean=-1;
%         data{1}.x=obj.T_0+1+obj.dy_T_mean*obj.z_list;
%         data{1}.y=obj.z_list;
%         plot_config.fontsize=20;
%         plot_config.label_list={1,'$\bar{T}_0+1+\bar{\mathcal{T}}_z z$', '$z$'};
%         plot_config.print_size=[1,500,900];
%         plot_config.name=['C:\Figure\DDC_LST\HB_porous_','T_0.png'];
%         plot_line(data,plot_config);
% 
%         data{1}.x=obj.d_T_0+obj.dy_T_mean;
%         data{1}.y=obj.z_list;
%         plot_config.fontsize=20;
%         plot_config.label_list={1,'$\partial_z \bar{T}_0+\bar{\mathcal{T}}_z$', '$z$'};
%         plot_config.print_size=[1,500,900];
%         plot_config.name=['C:\Figure\DDC_LST\HB_porous_','d_T_0.png'];
%         plot_line(data,plot_config);
%         
%         data{1}.x=obj.S_0+1+obj.dy_S_mean*obj.z_list;
%         data{1}.y=obj.z_list;
%         plot_config.label_list={1,'$\bar{S}_0+1+\bar{\mathcal{S}}_z z$', '$z$'};
%         plot_config.print_size=[1,500,900];
%         plot_config.name=['C:\Figure\DDC_LST\HB_porous_','S_0.png'];
%         plot_line(data,plot_config);
% 
%         data{1}.x=obj.S_0+1+obj.dy_S_mean*obj.z_list;
%         data{1}.y=obj.z_list;
%         plot_config.label_list={1,'$\bar{S}_0+1+\bar{\mathcal{S}}_z z$', '$z$'};
%         plot_config.print_size=[1,500,900];
%         plot_config.ylim_list=[1,0,0.01];
%         plot_config.name=['C:\Figure\DDC_LST\HB_porous_','S_0_local.png'];
%         plot_line(data,plot_config);
%         plot_config.ylim_list=0;
% 
%         data{1}.x=obj.S_0+1+obj.dy_S_mean*obj.z_list;
%         data{1}.y=obj.z_list;
%         plot_config.label_list={1,'$\bar{S}_0+1+\bar{\mathcal{S}}_z z$', '$z$'};
%         plot_config.print_size=[1,500,900];
%         plot_config.ylim_list=[1,0.49,0.51];
%         plot_config.xlim_list=[1,0.49,0.51];
%         plot_config.name=['C:\Figure\DDC_LST\HB_porous_','S_0_local_core.png'];
%         plot_line(data,plot_config);
%         plot_config.ylim_list=0; plot_config.xlim_list=0;
% 
% 
%         data{1}.x=obj.T_hat;
%         data{1}.y=obj.z_list;
%         plot_config.label_list={1,'$\widehat{T}$', '$z$'};
%         plot_config.print_size=[1,500,900];
%         plot_config.name=['C:\Figure\DDC_LST\HB_porous_','T_hat.png'];
%         plot_line(data,plot_config);
% 
%         data{1}.x=obj.S_hat;
%         data{1}.y=obj.z_list;
%         plot_config.label_list={1,'$\widehat{S}$', '$z$'};
%         plot_config.print_size=[1,500,900];
%         plot_config.name=['C:\Figure\DDC_LST\HB_porous_','S_hat.png'];
%         plot_line(data,plot_config);
% 
% 
%         data{1}.x=obj.S_hat;
%         data{1}.y=obj.z_list;
%         plot_config.label_list={1,'$\widehat{S}$', '$z$'};
%         plot_config.print_size=[1,500,900];
%         plot_config.ylim_list=[1,0,0.01];
%         plot_config.name=['C:\Figure\DDC_LST\HB_porous_','S_hat_local.png'];
%         plot_line(data,plot_config);
%         plot_config.ylim_list=0;
% 
%         data{1}.x=obj.w_hat;
%         data{1}.y=obj.z_list;
%         plot_config.label_list={1,'$\widehat{w}$', '$z$'};
%         plot_config.print_size=[1,500,900];
%         plot_config.name=['C:\Figure\DDC_LST\HB_porous_','w_hat.png'];
%         plot_line(data,plot_config);
% %     end
% end


