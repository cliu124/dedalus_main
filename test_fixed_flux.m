clear all;
close all;
clc;


kx=2*pi;
Ra_T_q_list_analytical=linspace(kx^4,50000,100);
w_max_analytical=2*sqrt((1-kx^4./Ra_T_q_list_analytical)./(2*kx^2./Ra_T_q_list_analytical));

Ra_T_q_list_DNS=[2000,4500,20000,30000,40000,50000];
w_max_DNS=[4.7291,12.207,30.566,37.959,44.13,49.539];

data{1}.x=Ra_T_q_list_analytical;
data{1}.y=w_max_analytical;
data{2}.x=Ra_T_q_list_DNS;
data{2}.y=w_max_DNS;
plot_config.label_list={1,'$Ra_{T,q}$','max($|w|$)'};
plot_config.Markerindex=3;
plot_config.user_color_style_marker_list={'k-','ro'};
plot_config.print_size=[1,1000,900];
plot_config.name='w_Ra_T_q_fixed_flux.png';
plot_line(data,plot_config);