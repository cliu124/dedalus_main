clear all;
close all;
clc;

Ra_T_q=[2*10^4;
        4*10^4;
        6*10^4;
        10^5;
        2*10^5;
        5*10^5;
        10^6;
        9*10^6];
    
    
Nu=[1.934513719434235;
    2.355410647995526;
    3.262885859488830;
    3.767677680127720;
    4.420324859766439;
    5.137144841275981;
    6.368201652567532;
    11.589667583344351];
Ra_T=Ra_T_q./Nu;
[eta,c0]=scaling(Ra_T(2:8),Nu(2:8));

clear all;
Ra_T_q=[10^8;
    3*10^8;
    10^9;
    3*10^9;
    10^10;
    ];
    
% Nu=[11.85;
%     12.03;
%     15.26;
%     22.67;
%     43.076
%     ];
Nu=[11.3095
12.8131
15.1509
20.1041
31.4387
];

Ra_T=Ra_T_q./Nu;
[eta_Ra_T,c0_Ra_T]=scaling(Ra_T,Nu);
[eta_Ra_T_q,c0_Ra_T_q]=scaling(Ra_T_q,Nu);

data{1}.x=Ra_T;
data{1}.y=Nu;
data{2}.x=Ra_T;
data{2}.y=c0_Ra_T*Ra_T.^eta_Ra_T;
plot_config.Markerindex=3;
plot_config.user_color_style_marker_list={'k*','k--'};
plot_config.name='Ra_T_Nu_scaling.png';
plot_config.loglog=[1,1];
plot_config.label_list={1,'$Ra_T$','$Nu$'};
plot_line(data,plot_config);

data{1}.x=Ra_T_q;
data{1}.y=Nu;
data{2}.x=Ra_T_q;
data{2}.y=c0_Ra_T_q*Ra_T_q.^eta_Ra_T_q;
plot_config.Markerindex=3;
plot_config.user_color_style_marker_list={'k*','k--'};
plot_config.loglog=[1,1];
plot_config.label_list={1,'$Ra_{T,q}$','$Nu$'};
plot_config.name='Ra_T_q_Nu_scaling.png';
plot_config.legend_list={1,'DNS','$Nu=0.189\;Ra_{T,q}^{0.217}$'};
plot_line(data,plot_config);
