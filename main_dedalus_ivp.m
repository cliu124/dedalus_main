clear all;
close all;
clc;

%%12073090: 32*32 size, elevator mode, without shear...
%%12073116: 96*32 size, elevator mode, without shear...
%%12075689: 96*32 size, elevator mode, with shear...
%%12075690: 32*32 size, elevator mode, with shear

%%These are local folder
% folder_name='C:\Data\dedalus\IFSC_2D_without_shear\';
%folder_name='C:\Data\dedalus\IFSC_2D_with_shear\';
% folder_name='C:\Data\dedalus\dedalus_12075689\IFSC_2D_with_shear\';
% folder_name='C:\Data\dedalus\dedalus_12080489\IFSC_2D_without_shear\';

% folder_name='C:\Data\dedalus\dedalus_12075689\IFSC_2D_with_shear\';

% folder_name='/rc_scratch/chli3324/dedalus_12073090/';
% file_name='IFSC_2D_without_shear_s1_random';
% file_name='IFSC_2D_without_shear_s1_small_domain';

% file_name='IFSC_2D_without_shear_s1_elevator_short';
% file_name='IFSC_2D_without_shear_s1_elevator_long';
% file_name='IFSC_2D_with_shear_s1';

% folder_name='C:\Data\dedalus\dedalus_12073090\IFSC_2D_without_shear\';
% file_name='IFSC_2D_without_shear_s1';
% folder_name='C:\Data\dedalus\dedalus_12075689\IFSC_2D_with_shear\';
% file_name='IFSC_2D_with_shear_s1';

slurm_num={'12073090',... %%IFSC, without shear, 32*32, A_elevator=1, A_noise=0
    '12073116',... %%IFSC, without shear, 32*96, A_elevator=1, A_noise=0
    '12075689',... %%IFSC, with shear, 32*96, A_elevator=1, A_noise=0.01
    '12083149',...%%IFSC without shear, 32*32, A_elevator=1, A_noise=0.01, A_shear=1, show the effect of shear
    '12083150',...%%IFSC, without shear, 32*32, A_elevator=0, A_noise=0.01, A_shear=1, show the effect of only shear... then decay
    '12083221',...%%IFSC, without shear, 32*32, A_elevator=1, A_noise=0.01, show effect of random noise...
    '12083491',...IFSC, with shear, 96*32, A_elevator=1, A_noise=0.01, ks=0.05
    '12083494',...IFSC with shear, 96*32, A_elevator=1, A_noise=0.01, ks=0.0264
    '12084941',...%%IFSC, without shear, 192*32, A_elevator=1,
    '12085402',...%%IFSC, with shear, 32*32 box size, Ra=1.1, very large initial condition in elevator and shear
    '12085400',...%%IFSC, with shear, 32*96, A_elevator=1, A_noise=0.01, ks: 0.008806549460544114
    '12085401',...%%IFSC, with shear, 32*96, A_elevator=1, A_noise=0.01, ks: 0.017613098921088227
    '12085887',...%%IFSC, without shear, 32*192, A_elevator=1, A_noise=0, A_shear=0
    '12086951',...%%IFSC, without shear, 96*32 domain, A_elevator=1, A_noise=0, A_shear=9
    '12088536',...%%IFSC, with shear, 32*32 domain, A_elevator=1, A_noise=0, A_shear=0
    '12088537',...%%IFSC, with shear, 32*32 domain, A_elevator=1, A_noise=0, A_shear=1
    '12088673',...%%IFSC, without shear, 8*192 domain, A_elevator=1, A_noise=0, A_shear=0
    '12089740',...%%IFSC, with shear, 32*32, A_elevator=1, A_noise=0, A_shear=1, long time up to 20000
    '12089741',...%%IFSC, with shear, 32*96, A_elevator=1, A_noise=0.01
    '12090288',...%%IFSC, with shear 8*8, A_elevator=1, A_noise=0.01, A_shear=0
    '12099609',...%%IFSC, with shear 8*24, ks=2*2pi/Lz, A_elevator=1, A_noise=0.01, u_L=1, F_sin: 0.001240885014416157
    '12099641',...%%IFSC, with shear 8*24, ks=2*2pi/Lz, A_elevator=1, A_noise=0.01, u_L=0.1, F_sin: 0.0001240885014416157
    '12099647',...%%IFSC, with shear 8*24, ks=2*2pi/Lz, A_elevator=1, A_noise=0.01, u_L=0.01, F_sin: 0.00001240885014416157
    '12099955',...%%IFSC, with shear 8*24, ks=2pi/Lz, A_elevator=1, A_noise=0.01, u_L=0.01, F_sin: 3.1022125360403927e-06, 
    '12099956',...%%IFSC with shear 8*24, ks=2pi/Lz, A_elevator=1, A_noise=0.01, u_L=0.1, F_sin: 3.1022125360403927e-05, 
    '12099957',...%%IFSC, with shear 8*24, ks=2pi/Lz, A_elevator=1, A_noise=0.01, u_L=1, F_sin: 3.1022125360403927e-04, 
    '12100768',...%%IFSC, with shear 8*24, ks=2pi/Lz, A_elevator=1, A_noise=0.01, u_L=1, F_sin: 3.1022125360403927e-04, F_sin_2ks: 0.001240885014416157, 
    '12100769',...%%IFSC, with shear 8*24, ks=2pi/Lz, A_elevator=1, A_noise=0.01, u_L=1, F_sin: 0.00031022125360403925, F_sin_2ks: 0.001240885014416157, F_sin_3ks: 0.0027919912824363536, 
    '12100793',...%%IFSC, with shear 8*24, ks=2pi/Lz, A_elevator=1, A_noise=0.01, u_L=1, F_sin: 0.00031022125360403925, F_sin_2ks: 0.001240885014416157, F_sin_3ks: 0.0027919912824363536, F_sin_4ks: 0.004963540057664628, 
    '12100794',...%%IFSC, with shear 8*24, ks=2pi/Lz, A_elevator=1, A_noise=0.01, u_L=1, F_sin: 3.1022125360403927e-04, F_sin_2ks: 0.001240885014416157, phase_2ks: pi/2
    '12100795',...%%IFSC, with shear 8*24, ks=2pi/Lz, A_elevator=1, A_noise=0.01, u_L=1, F_sin: 0.00031022125360403925, F_sin_2ks: 0.001240885014416157, phase_2ks: pi/4
    '12131253',...%%IFSC, with shear, 8*24, ks: 0.4227143741061175, A_elevator=1, A_noise=0, u_L=1,
    '12131254',...%%IFSC, with shear, 8*24, ks: 0.10567859352652938, A_elevator=1, A_noise=0, u_L=1,
    '12131256',...%%IFSC, with shear, 8*24, ks: 0.07045239568435291, A_elevator=1, A_noise=0, u_L=1,
    '12131257',...%%IFSC, with shear, 8*24, ks: 0.21135718705305875, A_elevator=1, A_noise=0, u_L=1,
    '12131258',...%%IFSC, with shear, 8*24, ks: 0.05283929676326469, A_elevator=1, A_noise=0, u_L=1,
    '12131259',...%%IFSC, with shear, 8*24, ks: 0.035226197842176454, A_elevator=1, A_noise=0, u_L=1,
    '12131260',...%%IFSC, with shear, 8*24, ks: 0.017613098921088227, A_elevator=1, A_noise=0, u_L=1,
    '12132188',...%%IFSC, with shear, 8*24, ks: 0.017613098921088227, A_elevator=1, A_noise=0.01, u_L=1,
    '12132197',...%%IFSC, with shear, 8*24, ks: 0.017613098921088227, A_elevator=1, A_noise=0.01, u_L=10^{-5}, F_sin: 3.1022125360403926e-09, 
    '12132206',...%%IFSC, with shear, 8*24, ks: 0.017613098921088227, A_elevator=1, A_noise=0.01, u_L=10^{-4}, F_sin: 3.1022125360403926e-08, 
    '12132211',...%%IFSC, with shear 8*24, ks: 0.017613098921088227,A_elevator=1, A_noise=0.01, u_L=10^{-3}, F_sin: 3.1022125360403926e-07, 
    '12132414',...%%IFSC, with shear 8*24, ks: 0.017613098921088227,A_elevator=1, A_noise=0.01, u_L=1, Ra_ratio=5
    '12132615',...%%IFSC, with shear, 8*24, ks: 0.017613098921088227,A_elevator=1, A_noise=0.01, u_L=10, Ra_ratio=1.1
    '12132761',...%%IFSC, with shear, 8*24,  ks: 0.017613098921088227,A_elevator=1, A_noise=0.01, u_L=10, Ra_ratio=2, time up to 1000
    '12132764',...%%IFSC, with shear, 8*24, ks: 0.017613098921088227,A_elevator=1, A_noise=0.01, u_L=1, Ra_ratio=5, time up to 500
    '12135159',...%%IFSC, with shear, 8*24, ks=0.0176, A_elevator=1, A_noise=0.01, u_
    '12135442',...%%IFSC, without shear, 8*24, Ra_ratio=5
    '12135952',...
    '12135952',...
    '12136034',...
    '12136695',...
    '12144003',...
    '12148188',...
    '12148389',...
    '12148590',...
    '12149388',...
    '12149389',...
    '12149390',...
    '12149391',...
    '12150307',...
    '12150308',...
    '12150309',...
    '12150310',...
    '12170693',...
    '12170695',...
    '12170697',...
    '12173787',...
    '12173788',...
    '12173789',...
    '12210887',...%%primitive, tau=0.01, Pr=10, R_rho_T2S=2, fingering, Lx=50, Lz=50
    '12211545',...%%primitive, tau=0.01, Pr=10, R_rho_T2S=0.99, diffusive, Lx=50, Lz=50
    '12213796',...%%primitive, tau=0.01, Pr=10, R_rho_T2S=0.99, diffusive, Lx=50, Lz=50
    '12213985',...%%primitive, tau=0, Pr=10, R_rho_T2S=0.99, diffusive, Lx=50, Lz=50, T=100
    '12213986',...%%primitive, tau=0, Pr=10, R_rho_T2S=2, fingering, Lx=50, Lz=50, T=100
    '12234192',...%%primitive, tau=0, Pr=10, R_rho_T2S=0.99, diffusive, Lx=50, Lz=50
    '12247548',...%%primitive, tau=0.01, Pr=10, R_rho_T2S=0.99, diffusive, Lx=50, Lz=50
    '12247549',...%%primitive, tau=0.01, Pr=10, R_rho_T2S=2, fingering, Lx=50, Lz=50
    '12248499',...%%primitive, tau=0.01, Pr=10, R_rho_T2S=0.5, diffusive, Lx: 1193.2458743605148, Lz: 18.644466786883044, with shear
    '12248567',...%%primitive, tau=0.01, Pr=10, R_rho_T2S=0.5, diffusive, Lx: 1193.2458743605148, Lz: 18.644466786883044, with shear, longer time
    '12263592',...
    '12269191',...
    '12269732',...
    '12288874',...
    '12288876',...
    '12288904',...
    '12288973',...
    '12288974',...
    '12289163',...
    '12289180',...
    '12289181',...
    '12291686',...
    '12291687',...
    '12291735',...
    '12291736',...
    '12291834',...
    '12292012',...
    '12292013',...
    '12356725',...
    '12356727',...'12357443',...'12357448',...%     '12357464',...
    '12363658',...
    'end'};
%     '12089742',...
%---Update simulation results for the bounded salt-finger to validate the
%staircase solutions
slurm_num={'13294812',... %For R_rho_T2S=40, tau=0.01, S1 IC
            '13294813',... %For R_rho_T2S=40, tau=0.01, S2 IC
            '13294814'... %For R_rho_T2S=40, tau=0.01, S3 IC
            };

slurm_num={'13294930',...
            '13294931'};
        
%%run 2022/03/16, put a range of initial conditions...
%%THese wavenumber setting is wrong... but can see some results there... 
% slurm_num={'13298996','13299000',...%Ra_S2T=1300,  %,'13299001'...
%         '13299003','13299004',... Ra_S2T=1400
%         '13299009','13299010',... %Ra_S2T=1500
%         '13299011','13299015',... %Ra_S2T=1600
%         '13299016','13299017',... %Ra_S2T=1700
%         '13299019','13299023',... %Ra_S2T=1800
%         '13299024','13299028','13299029', ...%Ra_S2T=1900
%         '13299054','13299056','13299062', ...%Ra_S2T=2000
%         '13299084','13299088','13299093', ...%Ra_S2T=2100
%         '13299136','13299137','13299138', ...%Ra_S2T=2200
%         '13299139','13299140','13299142', ...%Ra_S2T=2300 
%         '13299143','13299144','13299145', ...%Ra_S2T=2400
%         '13299150','13299151','13299152',... %Ra_S2T=2500
%         '13299153','13299154','13299155',...%Ra_S2T=2600
%         };

slurm_num={'13298996','13299000',...%Ra_S2T=1300,  %,'13299001'...
        '13299003','13299004',... Ra_S2T=1400
        '13299009','13299010',... %Ra_S2T=1500
        '13299011','13299015',... %Ra_S2T=1600
        '13299016','13299017',... %Ra_S2T=1700
        '13299019','13299023',... %Ra_S2T=1800
            };
slurm_num={'13299009','13299010',... Ra_S2T=1500
        '13299054','13299056','13299062'...Ra_S2T=2000
        };
%The correct wavenumber pair and domain size Lx2d=4, four finger pairs
% slurm_num={'13299208','13299209'}
slurm_num={'13299208','13299209',... Ra_S2T=1500
        '13299210','13299211','13299212',...Ra_S2T=2000
        '13299322','13299323','13299324',...Ra_S2T=2500
        '13299325','13299326','13299327',...Ra_S2T=3000
        };
    
% 
%  %The correct wavenumber pair and domain size Lx2d=1, one finger pair
% slurm_num={'13299383','13299384',... Ra_S2T=1500
%            '13299385','13299387','13299388',... Ra_S2T=2000
%            '13299389','13299390','13299391',... Ra_S2T=2500
%            '13299393','13299394','13299395',... Ra_S2T=3000
%            };

slurm_num={'13299932'};

%check the stability of periodic solutions.
% slurm_num={'13379372','13379373','13379374'};
slurm_num={'13379264','13379264'};
% slurm_num={'13379734','13379735',...
%     '13379740','13379741','13379750'};
% slurm_num={'13379742'};

%for the domain size Lx2d=1, and then for each Ra, put S1, S2, S3 solution
%profile as the initial conditions...
slurm_num={'13410044','13410045',...% '13410046', ...Ra_S2T=1/90*Ra_T
    '13410047','13410049','13410050',...Ra_S2T=1/80*Ra_T
    '13410052','13410053','13410054',...Ra_S2T=1/70*Ra_T
    '13413728','13413730','13413731',......Ra_S2T=1/60*Ra_T
    '13410299','13410300','13410303',...Ra_S2T=1/50*Ra_T
    '13398772', '13398773','13398774',...Ra_S2T=2500, Ra_S2T=1/40*Ra_T
    '13410250','13410251','13410256',...Ra_S2T=1/30*Ra_T,
    '13410263','13410264','13410265',...Ra_S2T=1/20*Ra_T
    '13398794','13398795','13404915'...Ra_S2T=10000
    '13431231'...Ra_S2T=20000
    %'13431232'...Ra_S2T=50000, need redo.. require small dt.
    %'13431234'...Ra_S2T=100000, need redo.. require small dt.
    };

%slurm_num=slurm_num(end);%This show evidence of traveling wave...

%tau=0.01, initialized by the S1 solution at R_rho_T2S=10...
% slurm_num={'13435676',...R_rho_T2S=19
%             '13435677',...R_rho_T2S=17
%             '13435678',...R_rho_T2S=15
%             '13435679',...R_rho_T2S=13
%             '13435680'...R_rho_T2S=11
%             };
% slurm_num=slurm_num(4);
%try to validate the nusselt nubmer against Yang's value
% slurm_num={'13435681'} %initialized by random noise, Lx=2.5
        
% slurm_num=slurm_num(end);

%tau=1/3, R_rho_T2S=1, initialized by the periodic solution from Hopf
%bifurcation of branch 2, and it is unstable to become the traveling wave.
% slurm_num={'13379264'}; %This is the simulation that show the periodic solution...
% slurm_num={'13431061'}; %run longer time... for traveling wave...

% slurm_num={'13415045'};%This is R_rho=10, tau=0.01, starting from S1 solution, Lx2d=1, with fine time sampling...
% slurm_num={'13435394'};This is R_rho=10, tau=0.01, start from 1 layer

%for the domain size Lx2d=4, check the nusselt number for 1-layer solution
% slurm_num={'13399059',... Ra_S2T=2500
%         '13399767', ... or '13399767'
%          '13399753' ...Ra_S2T=20000
%          '13399765'...Ra_S2T=50000
%          '13399766'...Ra_S2T=100000
%         }
    
    
% slurm_num=slurm_num(3:end);
%         ,...Ra_S2T=10000
%         '13399072', or
%         '13398982', or 
%         '13399100', or 
%         ''};

%2022/04/04, tau=0.3, test the initial condition as steady tilted convection roll
%The first column corresponds to the initial condition as steady symmetric
%roll. The second column is the steady tilted roll. 
slurm_num={'13447622',...R_rho_T2S=0.5
            '13448773','13447623',...R_rho_T2S=0.6
            '13447624',...R_rho_T2S=0.7
            '13447812',...R_rho_T2S=0.8
            '13447621',...R_rho_T2S=0.9
            '13447817'...R_rho_T2S=1
            }
        
%%2022/04/07, the initial condition as the stable tilted convection roll
%for tau=0.01, Pr=7 case
slurm_num={'13463098',...Ra_S2T=3258
           '13463112',...Ra_S2T=5542
           '13463117',...Ra_S2T=8766
           '13463118',...Ra_S2T=16290
           '13463119',...Ra_S2T=33487};
           '13463122'...Ra_S2T=59281
           };
 slurm_num=slurm_num(6);

%%2022/04/28 tau=0.01, R_rho_S2T=40, Ra_S2T=2500
slurm_num={'13562813', ...%Pr=7
        '13562812'... Pr=0.03
        };
slurm_num=slurm_num(2);

%%2022/05/01, wavenumber continuation
slurm_num={'13619823',... kx=19
            '13619824',... kx=18
            '13619825',... kx=16
            '13619826',... kx=15
            '13619827',... %kx=13
            '13619852', ...%kx=11
            '13619853',...%kx=9
            '13619854',...%kx=7
            '13619855',...%kx=5
            '13619856',...%kx=3
            '13619857',...%kx=1
            };
% slurm_num=slurm_num(end-1:end);
%2022/05/02 results for even wavenumber
%Pr=7, tau=0.01, Ra_T=10^5, R_rho=2, bounded salt finger
slurm_num={'13623956',...: kx=1
            '13623957',...: kx=2
            '13623959',...: kx=4
            '13623960',...: kx=6
            '13623961',...: kx=8
            '13623962',...: kx=10
            '13623963',...: kx=12
            '13623964',...: kx=14
            '13623965',...: kx=16
            '13623966'...: kx=18
            };
% 
% %even wavenumber, Pr=0.05, initial condition as S1
slurm_num={'13623969',...: kx=1
            '13623970',...: kx=2
            '13623971',...: kx=4
            '13623972',...: kx=6
            '13623973',...: kx=8
            '13623974',...: kx=10
            '13623975',...: kx=12
            '13623976',...: kx=14
            '13623977',...: kx=16
            '13623978'...: kx=18
            };
%Pr=0.05, initial condition as Tilted finger 1.
slurm_num={'13625602',...: kx=2
        '13625603',...: kx=4
        '13625604',...: kx=6
        '13625606',...: kx=8
        '13625607',...: kx=10
        '13625643',...: kx=12
        '13626014',...: kx=14
        '13626015'...: kx=16
        };
    
%Pr=0.05, initial condition as the traveling wave
% slurm_num={'13626016',...: kx=3
%             '13626034',...: kx=6
%             '13626035',...: kx=9
%             '13626037',...: kx=12
%             '13626131'...: kx=15
%             };
slurm_num={'13627393',...
            '13627394'};
slurm_num={'13629367'};
slurm_num={'13633501',...
            '13633505',...
            '13633506'};
slurm_num=slurm_num(2:3);
slurm_num={'13623961'};

%Pr=7, with IC from S2
slurm_num={'13635931',... kx=2
        '13635933',... kx=4
        '13635934',... kx=6
        '13635936',... kx=8:
        '13635937',... kx=10:
        '13635939',... kx=12:
        '13635940',... kx=14:
        '13635941',... kx=16:
        };
slurm_num={'13635936'};

%IC: S3
slurm_num={'13635943',... kx=4:
            '13635944',... kx=6:
            '13635946',... kx=8:
            '13635947',... kx=10:
            '13635949',... kx=12:
            '13635950',... kx=14 
            };

%Pr=0.05, TW1
slurm_num={'13635953',... kx=4:
            '13635954',... kx=6:
            '13635956',... kx=8:
            '13635957',... kx=10:
            '13635958',... kx=12:
            '13635959',... kx=14:
            '13635960'... kx=16:
            };
% slurm_num=slurm_num(3);
% slurm_num={'13635936'};
%IC from TF1, Pr=7, R_rho_T2S=40
% slurm_num={'13639120',... kx=4
%             '13639121'... kx=6
%             };
% slurm_num={'13639230'};
slurm_num={'13640186'};
slurm_num={'13908442','13908443','13908445','13908449'};
slurm_num={'13908443'};
slurm_num={'13910505',...: kx=18
            '13910508',...: kx=16
            '13910509',...: kx=14
            '13910874'...kx=12
            };
slurm_num={'13910922',...kx=10
    '13910962',... kx=8
    '13911003',... kx=6
    '13911062'...kx=4
    };
slurm_num={'13912354'};
        %,...: kx=10'};
% slurm_num={'13910236'};
slurm_num={'14170492',...
            '14170493',...
            '14170494',...
            '14170495',...
            '14170496',...
            '14170497',...
            '14170498',...
            '14170499',...
            '14170500'};
        
% slurm_num={'14172733',...
%             '14172834',...
%             '14173006',...
%             '14173419',...
%             '14174697'};
% slurm_num={'14175141'};
slurm_num={'14175555'};
% slurm_num={'14319541'};
slurm_num={'14319532',...
            '14319533',...
            '14319536',...
            '14319541',...
            '14319543',...
            '14319544',...
            '14319546'};
slurm_num={'14380195',...
            '14380196',...
            '14380197',...
            '14380198',...
            '14380209',...
            '14380210'};
% slurm_num={'14380197'};
% slurm_num={'14380210'};
slurm_num={'14380404',...
            '14380405',...
            '14380406',...
            '14380407'};
slurm_num={'14447605',...
            '14447606',...
            '14447607',...
            '14447608',...
            '14447609',...
            '14447611'};
slurm_num={'14450461',...
           '14450462',...
           '14450463',...
           '14450464',...
           '14450466',...
           '14450467'};
% slurm_num={'14450467'};
slurm_num={'14508196',...
            '14508746',...
            '14508750',...
            '14508752',...
            '14508753'};
% slurm_num={'14541125',...
%             '14541126'};
slurm_num={'14541127',...
            '14541128'};
% slurm_num={'14591433'};
slurm_num={'13635946'};
slurm_num={'14671572',...
    '14671573',...
    '14671574'};
% slurm_num={'14671567'};
slurm_num={'14673046',...
            '14673047',...
            '14673048',...
            '14673049'};
slurm_num={'14673047'};
slurm_num={'14672941',...
        '14672942',...
        '14672943'};
slurm_num={'14680804',...
    '14680806',...
    '14680814',...
    '14680810',...
    '14680791'};
slurm_num={'14739403',...
    '14739404',...
    '14739405',...
    '14739406',...
    '14739407'};
slurm_num={'14739406'};
        %,...
         %   '14672941',...
         %   '14672942',...
         %   '14672943'};
% slurm_num={'14672941',...
%             '14672942',...
%             '14672943'};
%         %            '14673894',...
slurm_num={'14741864',...
            '14741865',...
            '14741866',...
            '14741867'};
slurm_num={'14758590'};
slurm_num={'14680791'};
slurm_num={'14741869'};
slurm_num={'14741864'};
slurm_num={'14768963','14768964'};
% slurm_num={'14741864','14741865'};
slurm_num={'14673046'}; %14739405
slurm_num={'14739406'};
slurm_num={'14739406'}; 
slurm_num={'15144692'};
% slurm_num={'15155150','15155156','15155162'};
% slurm_num={'15192717'};
% slurm_num={}
%14768964, 14741865, 14768963
slurm_num={'907867',
'907878',
'15154914',
'15156185',
'15156187',
'15156193',
'15156832',
'15156833',
'15156834',
'15156835',
'15157523',
'15159128',
'15159132',
'15159726',
'15163269',
'15163270'
};
slurm_num={'15155309','15155314'};
slurm_num={'15295430'};
slurm_num={'15313466'};
slurm_num={'15373305',...
    '15373308'};

slurm_num={'15376313',
'15376314',
'15376315',
'15373305',
'15373308',
'15376324',
'15376325',
'15376327',
'15376328',
'15376323'
};

slurm_num={'15380122',
'15380125',
'15380128',
'15380131',
'15380134',
'15380137',
'15380138'
};

slurm_num={'15245854',
'15287990',
'15287991',
'15287998',
'14739406',
'15295286',
'907867',
'907878',
'15154914',
'15156185',
'15156187',
'15156193',
'15156832',
'15156833',
'15156834',
'15156835',
'15157523',
'15159128',
'15159132',
'15159726',
'15163269',
'15163270',
'15263842',
'15163950',
'15163951',
'15163952',
'15163953',
'15163271',
'15159728',
'15157545',
'15156198',
'15154920',
'15154925',
'907903',
'907919',
'907934',
'907946',
'15295436',
'15295437',
'14739405'
    };
slurm_num={'15394862',
'15376324',
'15376325',
'15376327',
'15376328',
'15376323',
'15379741',
'15379742',
'15379747',
'15379750',
'15379754'
};
slurm_num={'15565226'};

slurm_num={'15597576',...
    '15597577',...
    '15597578',...
    '15597579'};
% slurm_num={};
% slurm_num=slurm_num(end);
flag.print=1; 
flag.video=0;
flag.visible=1;
flag.no_ylabel=0;
% 
% 
% for slurm_ind=1:length(slurm_num)
%     h5_name=['D:\Data\dedalus\dedalus_',...
%         slurm_num{slurm_ind},...
%         '\analysis\analysis_s1.h5'];
% 
%      set(0,'DefaultFigureVisible','on')
%      dedalus_post_my{slurm_ind}=dedalus_post(h5_name,flag);
%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.dedalus_post_ivp();
%      dedalus_post_my{slurm_ind}.print=0; dedalus_post_my{slurm_ind}.visible=0;
%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.spectrum_t('u_x_ave',[0.25],[],[2]);
%      
%      %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('u');
%      data{1}.y(slurm_ind)=dedalus_post_my{slurm_ind}.freq_sort(1);
%      data{1}.x(slurm_ind)=dedalus_post_my{slurm_ind}.Ra_T;
%      
% end
% data{1}.y=data{1}.y*2*pi;
% plot_config.label_list={1,'$Ra_{T,q}$','$\omega$'};
% plot_config.name='RBC_Ra_global_Ra_Tq_omega.png';
% plot_line(data,plot_config);

% data{2}.x=data{1}.x;
% data{2}.y=1./(-log(abs(data{2}.x-46892.05)));

% error('1');
for slurm_ind=1:length(slurm_num)
    h5_name=['D:\Data\dedalus\dedalus_',...
        slurm_num{slurm_ind},...
        '\analysis\analysis_s1.h5'];

     set(0,'DefaultFigureVisible','on')
     dedalus_post_my{slurm_ind}=dedalus_post(h5_name,flag);
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('dy_T_mean_q');   

     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.z_slice('T',[0.1,0.3,0.5],[1,300]);

     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('u');   


     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.total_xt_ave('T',[],[276,466]);


     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.get_Nu('T',[1,1000]);

     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.dedalus_post_ivp();
     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('u');
     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('T');
     
     % error('1');
     %dedalus_post_my{slurm_ind}.title_time=0;
     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.snapshot('T',[3],[],[1000]);
%      error('1');
%      error('1');
     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.spectrum_snapshot('T');
     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.spectrum_average('T',[0.1,0.3]);
%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.total_xt_ave('u',[],[276,466]);

     %      error('1');

%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.get_Nu('T',[500]);
%     error('1');
%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('u',[]);
%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('T');
%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('u');
%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.spectrum_t('Nu_T_t',[0.3],[],[100,500]);
%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.spectrum_t('u',[0.4],[],[]); %u

     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.phase_diagram('u','T',[],'max_z2')
%      error('1');

     % 
%      dedalus_post_my{slurm_ind}.print=1; dedalus_post_my{slurm_ind}.visible=0;
%      dedalus_post_my{slurm_ind}.video=1;
%      dedalus_post_my{slurm_ind}.title_time=1;
end
error('1');
for slurm_ind=1:2
    data{slurm_ind}.x=dedalus_post_my{slurm_ind}.t_list;
    data{slurm_ind}.y=dedalus_post_my{slurm_ind}.Nu_T_t_full;
end
plot_config.Markerindex=3;
plot_config.user_color_style_marker_list={'k-','b--'};%,'r-.','m:'};
plot_config.label_list={1,'$t$','$nu(t)$'};
plot_config.legend_list={1,'$\beta=1$','$\beta=10$'};%,'$\beta=100$','$\beta=1000$'};
plot_config.name='Nu_t_Ra_Tq_6e4_beta_1_10.png';
plot_line(data,plot_config);

for slurm_ind=3:4
    data{slurm_ind-2}.x=dedalus_post_my{slurm_ind}.t_list;
    data{slurm_ind-2}.y=dedalus_post_my{slurm_ind}.Nu_T_t_full;
end
plot_config.Markerindex=3;
plot_config.user_color_style_marker_list={'k-','b--'};%,'r-.','m:'};
plot_config.label_list={1,'$t$','$nu(t)$'};
plot_config.legend_list={1,'$\beta=100$','$\beta=1000$'};%,'$\beta=100$','$\beta=1000$'};
plot_config.name='Nu_t_Ra_Tq_6e4_beta_100_1000.png';
plot_line(data,plot_config);



error('1');
for slurm_ind=1:length(slurm_num)
    Nu(slurm_ind,1)=dedalus_post_my{slurm_ind}.Nu;
end
for slurm_ind=1:length(slurm_num)%:length(slurm_num)-1%[find(strcmp(slurm_num,'12247549'))]%slurm_ind=length(slurm_num)-2:length(slurm_num)-1
    %find(strcmp(slurm_num,'12136034'))
    %length(slurm_num)-1:length(slurm_num)-1
    
    %%change the path into D... just store data in the external disk...
    h5_name=['D:\Data\dedalus\dedalus_',...
        slurm_num{slurm_ind},...
        '\analysis\analysis_s1.h5'];

     set(0,'DefaultFigureVisible','on')
     dedalus_post_my{slurm_ind}=dedalus_post(h5_name,flag);

     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.dedalus_post_ivp();
     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.E_time('T',0);
%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.E_time('T');

     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('T',[201,1000]);
     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('p');
     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('w');
     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('T');

     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('w');
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('w');
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('u');
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('ww');
     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave_z_max('u',[],[]);

     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('dy_T_mean_q');
     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.get_Nu('T',[150]);

%      dedalus_post_my{slurm_ind}.z_basis_mode='Chebyshev';
%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.spectrum_average('S');
%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.spectrum_average('u');
%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.spectrum_average('w');

%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('S');

%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('rho');

     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.u_fluctuation_x_ave();
%     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.spectrum_t('S',0.5,0.5,[20,30]);
%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.rms_xt('u');
%     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.spectrum_average('S');
%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.total_xt_ave('S');
%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.total_xt_ave('S',[],[]);
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.total_xt_ave('T',[],[21,261]);

     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.total_xt_ave('w');
%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.total_xt_ave('w');

     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.total_xt_ave('rho');

     dedalus_post_my{slurm_ind}.print=1; dedalus_post_my{slurm_ind}.visible=1;
     dedalus_post_my{slurm_ind}.video=0;
     dedalus_post_my{slurm_ind}.title_time=1;
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.snapshot('T',5);

     dedalus_post_my{slurm_ind}.print=1; dedalus_post_my{slurm_ind}.visible=0;
     dedalus_post_my{slurm_ind}.video=0;
     dedalus_post_my{slurm_ind}.title_time=0;
     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.snapshot('T',[],[],[385,391,399]);
%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.snapshot('T',[],[],[309,313,321]);
     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.snapshot('T',[],[],[126,161]);
     
     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.snapshot('T',[],[],590:610);

%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.snapshot('T',[],[],[1000]);

     
     dedalus_post_my{slurm_ind}.video=0;
     dedalus_post_my{slurm_ind}.print=0;
     dedalus_post_my{slurm_ind}.title_time=1;
     dedalus_post_my{slurm_ind}.no_ylabel=1;
%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.total_xt_ave('S',[],[80,100]);

     
     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.snapshot('w');

%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.snapshot('T');

     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.u_laminar();

end

switch slurm_num{1}
    case {'13640186','13639526','13639527','13742323'}
        clear data;
        switch slurm_num{1}
            case {'13640186'}
                %kx=10
                p_TF1=load(['C:\Data\pde2path\HB_benard_full_zonal\',...
                    'salt_finger_kx_low_Ra_S2T_low_Pr_2D_no_slip\tr\bpt1\bpt1\pt177.mat']);
                p_TW1=load(['C:\Data\pde2path\HB_benard_full_zonal\',...
                    'salt_finger_kx_low_Ra_S2T_low_Pr_2D_no_slip\tr\bpt1\bpt2\pt161.mat']);
                TW1_DNS_range=[81:181];
                TF1_DNS_range=[281:381];
                TF1_LSS_snapshot=[333,334,335];
            case '13639526'
                %kx=14;
                p_TF1=load(['C:\Data\pde2path\HB_benard_full_zonal\',...
                    'salt_finger_kx_low_Ra_S2T_low_Pr_2D_no_slip\tr\bpt1\bpt1\pt123.mat']);
                p_TW1=load(['C:\Data\pde2path\HB_benard_full_zonal\',...
                    'salt_finger_kx_low_Ra_S2T_low_Pr_2D_no_slip\tr\bpt1\bpt2\pt94.mat']);
                
                TW1_DNS_range=[101:301];
                TF1_DNS_range=[600:1000];
                TF1_LSS_snapshot=[600,601,602];
            case '13639527'
                %kx=16
                p_TF1=load(['C:\Data\pde2path\HB_benard_full_zonal\',...
                    'salt_finger_kx_low_Ra_S2T_low_Pr_2D_no_slip\tr\bpt1\bpt1\pt58.mat']);
                p_TW1=load(['C:\Data\pde2path\HB_benard_full_zonal\',...
                    'salt_finger_kx_low_Ra_S2T_low_Pr_2D_no_slip\tr\bpt1\bpt2\pt7.mat']);
                
                TW1_DNS_range=[80:180];
                TF1_DNS_range=[300:400];
                TF1_LSS_snapshot=[333,334,335];
            case '13742323'
                p_TF1=load(['C:\Data\pde2path\HB_benard_full_zonal\',...
                    'salt_finger_kx_low_Ra_S2T_2D_no_slip\tr\bpt1\pt217.mat']);
                p_TW1=load(['C:\Data\pde2path\HB_benard_full_zonal\',...
                    'salt_finger_kx_low_Ra_S2T_2D_no_slip\tr\bpt1\pt217.mat']);
                TW1_DNS_range=[60:80];
                TF1_DNS_range=[60:80];
                TF1_LSS_snapshot=[60,61,62];
                
        end
        
        data{1}.x=mean(mean(dedalus_post_my{1}.S(:,:,TW1_DNS_range),2),3)+dedalus_post_my{1}.z_list;
        data{1}.y=dedalus_post_my{1}.z_list;
        data{2}.x=mean(mean(dedalus_post_my{1}.S(:,:,TF1_DNS_range),2),3)+dedalus_post_my{1}.z_list;
        data{2}.y=dedalus_post_my{1}.z_list;
        data{3}.x=dedalus_post_my{1}.z_list;
        data{3}.y=dedalus_post_my{1}.z_list;
        plot_config.label_list={1,'$z+\langle S \rangle_{h,t}$','$z$'};
        plot_config.print_size=[1,500,900];
        plot_config.Markerindex=3;
        plot_config.user_color_style_marker_list={'k-','r--','b-.'};%cl_list(A2_before_ind+1:A2_before_ind+2);
%         plot_config.name=[my.folder_name,'profile_together_TF1_TW1','_par=',num2str(point_list(point_ind)),'_profile_',var_name,'.png'];
        plot_config.legend_list={1,['$t\in [',num2str(min(TW1_DNS_range)-1),',',num2str(max(TW1_DNS_range)-1),']$'],...
            ['$t\in [',num2str(min(TF1_DNS_range)-1),',',num2str(max(TF1_DNS_range)-1),']$'],'$z$'};
        plot_config.fontsize_legend=16;
        plot_config.fontsize=34;
        plot_config.xlim_list=0; plot_config.xtick_list=0;
        plot_config.linewidth=3;
        plot_config.name=[dedalus_post_my{1}.h5_name(1:end-3),'_','S_0','_total_xt_ave_OTF1_TW1.png'];
        plot_config.print_size=[1,500,900];
        plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
        plot_line(data,plot_config);

        %data{1}.x=mean(mean(dedalus_post_my{1}.S(:,:,80:180),2),3)+dedalus_post_my{1}.z_list;
        %data{1}.y=dedalus_post_my{1}.z_list;
        %add the comparison with single mode solution 2022/05/23
        
        clear data;
        data{1}.x=mean(mean(dedalus_post_my{1}.S(:,:,TF1_DNS_range),2),3)+dedalus_post_my{1}.z_list;
        data{1}.y=dedalus_post_my{1}.z_list;
        data{2}.x=p_TF1.p.u(1+4*p_TF1.p.np:5*p_TF1.p.np)+p_TF1.p.x;
        data{2}.y=p_TF1.p.x;
        data{3}.x=dedalus_post_my{1}.z_list;
        data{3}.y=dedalus_post_my{1}.z_list;
        plot_config.label_list={1,'$z+\langle S \rangle_{h,t}$','$z$'};
        plot_config.user_color_style_marker_list={'k-','r--','b-.'};%cl_list(A2_before_ind+1:A2_before_ind+2);
        plot_config.legend_list={1,['$t\in [',num2str(min(TF1_DNS_range)-1),',',num2str(max(TF1_DNS_range)-1),']$'],'$z+\bar{S}_0$ in SME', '$z$'};
        plot_config.name=[dedalus_post_my{1}.h5_name(1:end-3),'_','S_0','_total_xt_ave_OTF1.png'];
        plot_line(data,plot_config);
        
       clear data
        data{1}.x=mean(dedalus_post_my{1}.u(:,:,TF1_LSS_snapshot(1)),2);
        data{1}.y=dedalus_post_my{1}.z_list;
        data{2}.x=mean(dedalus_post_my{1}.u(:,:,TF1_LSS_snapshot(2)),2);
        data{2}.y=dedalus_post_my{1}.z_list;
        data{3}.x=mean(dedalus_post_my{1}.u(:,:,TF1_LSS_snapshot(3)),2);
        data{3}.y=dedalus_post_my{1}.z_list;
        data{4}.x=p_TF1.p.u(1+8*p_TF1.p.np:9*p_TF1.p.np);
        data{4}.y=p_TF1.p.x;
        plot_config.label_list={1,'$ \langle u \rangle_h$','$z$'};
        plot_config.user_color_style_marker_list={'k-','r--','m:','b-.'};%cl_list(A2_before_ind+1:A2_before_ind+2);
        plot_config.legend_list={1,['$t=',num2str(TF1_LSS_snapshot(1)-1),'$'],...
            ['$t=',num2str(TF1_LSS_snapshot(2)-1),'$'],...
            ['$t=',num2str(TF1_LSS_snapshot(3)-1),'$'],...
            '$\bar{U}_0$ in SME'};        
        plot_config.name=[dedalus_post_my{1}.h5_name(1:end-3),'_','U_0','_total_xt_ave_TF1.png'];
        plot_config.xlim_list=[0,-0.15,0.15];
        plot_line(data,plot_config);
        
        clear data;
        data{1}.x=mean(mean(dedalus_post_my{1}.S(:,:,TW1_DNS_range),2),3)+dedalus_post_my{1}.z_list;
        data{1}.y=dedalus_post_my{1}.z_list;
        data{2}.x=p_TW1.p.u(1+4*p_TW1.p.np:5*p_TW1.p.np)+p_TW1.p.x;
        data{2}.y=p_TW1.p.x;
        data{3}.x=dedalus_post_my{1}.z_list;
        data{3}.y=dedalus_post_my{1}.z_list;
        plot_config.label_list={1,'$ z+\langle S \rangle_{h,t}$','$z$'};
        plot_config.user_color_style_marker_list={'k-','r--','b-.'};%cl_list(A2_before_ind+1:A2_before_ind+2);
        plot_config.legend_list={1,['$t\in [',num2str(min(TW1_DNS_range)-1),',',num2str(max(TW1_DNS_range)-1),']$'],...
            '$z+\bar{S}_0$ in SME','$z$'};        
        plot_config.name=[dedalus_post_my{1}.h5_name(1:end-3),'_','S_0','_total_xt_ave_TW1.png'];
        plot_line(data,plot_config);
        
        
        clear data
        data{1}.x=mean(mean(dedalus_post_my{1}.u(:,:,TW1_DNS_range),2),3);
        data{1}.y=dedalus_post_my{1}.z_list;
        data{2}.x=p_TW1.p.u(1+8*p_TW1.p.np:9*p_TW1.p.np);
        data{2}.y=p_TW1.p.x;
        plot_config.label_list={1,'$\langle u \rangle_{h,t}$','$z$'};
        plot_config.user_color_style_marker_list={'k-','r--','b-.'};%cl_list(A2_before_ind+1:A2_before_ind+2);
        plot_config.legend_list={1,['$t\in [',num2str(min(TW1_DNS_range)-1),',',num2str(max(TW1_DNS_range)-1),']$'],'$\bar{U}_0$ in SME'};        
        plot_config.name=[dedalus_post_my{1}.h5_name(1:end-3),'_','U_0','_total_xt_ave_TW1.png'];
        plot_line(data,plot_config);
        
end

error('1')
% 
% vslow=VideoWriter('TravelingWaveTwo20181111.avi');
% vslow.FrameRate=10;
% open(vslow);
% writeVideo(vslow,Fxnorm);
% close(vslow);

    % file_name='';
    % h5_name=[folder_name,file_name,'.h5'];


    % u=h5read(h5_name,'/tasks/u');
    % w=h5read(h5_name,'/tasks/w');
    % T=h5read(h5_name,'/tasks/T');
%     S=h5read(h5_name,'/tasks/S');
%     S_coeff=h5read(h5_name,'/tasks/S_coeff');
%     S_coeff=S_coeff.r+1i*S_coeff.i;


Ra_ratio=1.1;
lambda_opt=sqrt(1/2*(-2-Ra_ratio+sqrt(Ra_ratio^2+8*Ra_ratio)))*(3*Ra_ratio-sqrt(Ra_ratio^2+8*Ra_ratio))/(sqrt(Ra_ratio^2+8*Ra_ratio)-Ra_ratio);
k_opt=(1/2*(-2-Ra_ratio+sqrt(Ra_ratio^2+8*Ra_ratio)))^(1/4);
ks=k_opt/32;
tau=0.01;
Pr=7;
R_rho=90.9;
Ri=1;
uL=1/tau/ks*sqrt(Pr*(1-1/R_rho)/Ri);


data{1}.x=kx_list/k_opt;
data{1}.y=kz_list(1:Nz/2)/k_opt;
data{1}.z=log10(mean(abs(S_coeff(1:Nz/2,:,:)),3));
plot_config.zlim_list=[1,-3,0];
plot_config.xlim_list=[1,0,2];
plot_config.ylim_list=[1,-2,2];
plot_config.ztick_list=[1,-3,-2,-1,0];
plot_config.print_size=[1,1100,900];
plot_config.label_list={1,'$k/k_{opt}$','$m/k_{opt}$'};
plot_config.colormap='white_zero';
plot_contour(data,plot_config);


Nx=length(x);
Nz=length(z);
for t_ind=1:length(t)
%     u_int(t_ind)=sum(sum(u(:,:,t_ind).^2));
    E.S(t_ind)=sum(sum(S(:,:,t_ind).^2))/Nx/Nz/2;
    %E.T(t_ind)=sum(sum(T(:,:,t_ind).^2))/Nx/Nz/2;
    %E.u(t_ind)=sum(sum(u(:,:,t_ind).^2))/Nx/Nz/2;
    %E.w(t_ind)=sum(sum(w(:,:,t_ind).^2))/Nx/Nz/2;
    %E.kinetic(t_ind)=E.u(t_ind)+E.w(t_ind);
    average.S(t_ind)=mean(mean(S(:,:,t_ind)));
    %average.T(t_ind)=mean(mean(T(:,:,t_ind)));
    %average.u(t_ind)=mean(mean(u(:,:,t_ind)));
    %average.w(t_ind)=mean(mean(w(:,:,t_ind)));

end

[val,max_ind]=max(E.S);
t_grow=t(1:max_ind);

data{1}.x=t;
data{1}.y=E.S;
% data{2}.x=t_grow;
% data{2}.y=E.S(max_ind)*exp(2*lambda_opt*(t_grow-max(t_grow)));
plot_config.label_list={1,'$t$','$E_S$'};
plot_config.legend_list={0};
% plot_config.legend_list={1,'Simulation','Linear stability'};
plot_config.name=[folder_name,file_name,'E_S.png'];
plot_config.Markerindex=3;
plot_config.user_color_style_marker_list={'k-','bo--'};
plot_line(data,plot_config);

% error('1')

% plot(t_grow,E.S(1)*exp(2*lambda_opt*t_grow)); hold on;
% plot(t,E.S);
% [~,t_ind_1200]=min(abs(t-1200));
clear data
