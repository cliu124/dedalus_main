clear all;
close all;
clc;

Ra_freq=[46759.00000000000000	11.1244
        46760.50000000000000	10.0968
        46760.80000000000000	9.5715
        46760.90000000000000	9.2886
        46761.05000000000000	8.3433
        46761.07500000000000	7.6032
        46761.08100000000000	6.8692
        46761.08175000000000	6.381
        46761.08190000000000	6.1022
        46761.08195000000000	5.7929
        46761.08197500000000	5.1308
        46761.08197575000000	4.9686
        46761.08197610000000	4.7603
        46761.08197620000000	4.5663
        46761.08197622000000	4.4688
        46761.08197623000000	4.3938
        46761.08197624000000	4.2
        46761.08197624250000	4.0148
        ];
data{1}.x=Ra_freq(:,1);
data{1}.y=2*pi./Ra_freq(:,2);
Ra_g=46761.08197624290000;
% Ra_g=46761.08197620000;
fit_ind=8:18;%:length(data{1}.x);
modelfun = @(b,x)(b(1).*(-log(Ra_g-x)));
mdl = fitnlm(data{1}.x(fit_ind),data{1}.y(fit_ind),modelfun,[0.07854]);
% error('1');
% x=lsqr(-log(Ra_g-data{1}.x(fit_ind)),data{1}.y(fit_ind));
data{2}.x=data{1}.x(fit_ind);
%old, that from fitting
data{2}.y=0.07854.*(-log(Ra_g-data{2}.x));

%new, try from theory
% Ra_g=46761.08197624290000;

% data{3}.x=data{1}.x;
% data{3}.y=(1/43.688).*(-log(Ra_g-data{1}.x))+1.0697;

% factor=data{1}.y(end-1)/data{2}.y(end-1);
% data{2}.y=factor*data{2}.y;
plot_config.legend_list={1,'DNS','$T_p=-0.07854\;{\rm ln}(Ra_{T,q}^{(g)}-Ra_{T,q})$'};
plot_config.label_list={1,'$Ra_{T,q}$','$T_p$'};
plot_config.name='RBC_Ra_global_SM_Ra_Tq_period.png';
plot_config.user_color_style_marker_list={'msquare','k--','b:'};
plot_config.Markerindex=3;
plot_config.print=1;
plot_config.linewidth=3;
plot_config.xlim_list=[1,Ra_g-0.01,Ra_g];
plot_line(data,plot_config);

clear plot_config data;
load('pt68.mat');
data{1}.x=-real(p.sol.muv);
data{1}.y=-imag(p.sol.muv);
% data{2}.x=0;
% data{2}.y=0;
plot_config.label_list={1,'Re($\lambda$)','Im($\lambda$)'};
% plot_config.label_list={1,'$Ra_{T,q}$','$T_p$'};
plot_config.name='eigenvalue_secondary_instability_elevator.png';
plot_config.user_color_style_marker_list={'k*'};
plot_config.Markerindex=3;
plot_config.print=1;
plot_config.xlim_list=[1,-150,50];
plot_config.linewidth=3;
plot_line(data,plot_config);

error('1');
Ra_g=46761.08197624300000;
data_T{1}.x=data{1}.x;
data_T{1}.y=2*pi./data{1}.y;
data_T{2}.x=data{2}.x;
data_T{2}.y=1/41.108.*(-log(Ra_g-data{2}.x));
plot_config.legend_list={1,'DNS','$T=0.024326[-{\rm ln}(Ra_{T,q}^{(g)}-Ra_{T,q})]$'};
plot_config.label_list={1,'$Ra_{T,q}$','$T$'};
plot_config.name='RBC_Ra_global_SM_Ra_Tq_T.png';
plot_config.user_color_style_marker_list={'msquare','k--'};
plot_config.Markerindex=3;
plot_config.print=1;
plot_config.linewidth=3;
plot_config.xlim_list=[1,Ra_g-0.01,Ra_g];
plot_line(data_T,plot_config);