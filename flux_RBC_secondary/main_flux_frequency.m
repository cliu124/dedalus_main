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
data{1}.y=Ra_freq(:,2);
Ra_g=46761.08197624290000;
fit_ind=8:length(data{1}.x);
modelfun = @(b,x)(b(1)./(-log(b(2)-x)));
mdl = fitnlm(data{1}.x(fit_ind),data{1}.y(fit_ind),modelfun,[80,Ra_g]);

data{2}.x=data{1}.x(fit_ind);
data{2}.y=80./(-log(Ra_g-data{2}.x));
% factor=data{1}.y(end-1)/data{2}.y(end-1);
% data{2}.y=factor*data{2}.y;
plot_config.legend_list={1,'DNS','$\omega=80/[-{\rm ln}(Ra_{T,q}^{(g)}-Ra_{T,q})]$'};
plot_config.label_list={1,'$Ra_{T,q}$','$\omega$'};
plot_config.name='RBC_Ra_global_SM_Ra_Tq_omega.png';
plot_config.user_color_style_marker_list={'msquare','k--'};
plot_config.Markerindex=3;
plot_config.print=1;
plot_config.linewidth=3;
plot_config.xlim_list=[1,Ra_g-0.01,Ra_g];
plot_line(data,plot_config);
