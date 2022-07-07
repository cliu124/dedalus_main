clear all;
close all;
clc;
table=[1  334.7 -1.282 12.82 -0.6668 -0.5421 4.973 -0.2985  243;
        2  167.36  -2.526 12.75 -0.0653 0.5588 8.028 0.2694  93;
        3 83.6 -1.694 6.545 -0.19212 -0.7571 3.140 -0.1999   82;
        4  33.47  -2.038 27.48  -0.3613 -0.05530 3.035 -0.3541 607;
        5  16.74 -1.403 6.453 -0.2150 -0.07924 2.030 0.3777  165;
        6  8.36  -1.983 10.94 -0.8262 -0.3124 2.931 -1.934 275;
        7  3.47  -3.057 4.790 -0.2161 2.607 7.689 -0.3367 164;
        8 1.67   -3.261 7.912 -0.2662 -0.05444 2.210 -0.078 95];
data{1}.x=table(:,2);
data{2}.x=table(:,2);
data{3}.x=table(:,2);
data{4}.x=table(:,2);
data{5}.x=table(:,2);
data{6}.x=table(:,2);

data{1}.y=table(:,3);
data{2}.y=table(:,4);
data{3}.y=table(:,5);
data{4}.y=table(:,6);
data{5}.y=table(:,7);
data{6}.y=table(:,8);

plot_config.label_list={1,'$R_\rho$',''};
plot_config.legend_list={1,'mean($w$)','var($w$)','skew($w$)',...
    'mean($u$)','var($u$)','skew($u$)'};
plot_config.user_color_style_marker_list={'k*-','ro-','b^-','k--*','r--o','b--^'};
plot_config.Markerindex=3;
plot_config.print_size=[1,900,900];
plot_config.fontsize_legend=18;
plot_config.name='u_w_R_rho.png';
plot_line(data,plot_config);

plot_config.xlim_list=[1,0,10];
plot_config.name='u_w_R_rho_small.png';
plot_line(data,plot_config);

