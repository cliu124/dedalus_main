clear all;
close all;
clc;
global Ra kx kxpp darcy_hewitt
darcy_hewitt='3D';
solver='two_mode';%{‘two_mode’,'one_mode'}
switch darcy_hewitt
    case '2D'
        Ra_list=[10000,20000,40000];
        kx_list=0.48*Ra_list.^0.4;
%         kx_list=0.5*Ra_list.^(0.5);
        plot_config.legend_list={1,'Ra=10000','Ra=20000','Ra=40000'};
        ratio_u_dw=1./kx_list;
    case '3D'
        Ra_list=[4000,8000,16000];
        kx_list=0.17*Ra_list.^0.52;
%        kxpp_list=[];
        kxpp_list=[79,137,207];
        plot_config.legend_list={1,'Ra=4000','Ra=8000','Ra=16000'};
        ratio_u_dw=kx_list./(2*kx_list.^2);
end
% Ra=40000;%[10000,20000,40000];
for Ra_ind=1:length(Ra_list)
    Ra=Ra_list(Ra_ind);
    kx=kx_list(Ra_ind);
    kxpp=kxpp_list(Ra_ind);
    xmesh=linspace(0,1,200);
    if strcmp(solver,'two_mode')
        solinit=bvpinit(xmesh,@guess);
        switch darcy_hewitt
            case '2D'
                sol=bvp4c(@bvpfcn_2D,@bcfcn,solinit);
            case '3D'
                sol=bvp4c(@bvpfcn_3D,@bcfcn,solinit);
        end
        solinit=bvpinit(sol.x,@guess_two_mode);
        solinit.y(1:6,:)=sol.y(1:6,:);
        switch darcy_hewitt
            case '2D'
                sol=bvp4c(@bvpfcn_2D_two_mode,@bcfcn_two_mode,solinit);
            case '3D'
                sol=bvp4c(@bvpfcn_3D_two_mode,@bcfcn_two_mode,solinit);
        end
    elseif strcmp(solver,'one_mode')
        solinit=bvpinit(xmesh,@guess);
        switch darcy_hewitt
            case '2D'
                sol=bvp4c(@bvpfcn_2D,@bcfcn,solinit);
            case '3D'
                sol=bvp4c(@bvpfcn_3D,@bcfcn,solinit);
        end
    end
    
    data_w_hat_rms{Ra_ind}.x=sol.x;
    data_w_hat_rms{Ra_ind}.y=sqrt(2)*sol.y(1,:);
    data_T_hat_rms{Ra_ind}.x=sol.x;
    data_T_hat_rms{Ra_ind}.y=sqrt(2)*sol.y(3,:);
    data_T0_bar{Ra_ind}.y=sol.x;
    data_T0_bar{Ra_ind}.x=sol.y(5,:);
    data_u_hat_rms{Ra_ind}.x=sol.x;
    data_u_hat_rms{Ra_ind}.y=sqrt(2)*abs(sol.y(2,:))*ratio_u_dw(Ra_ind);
    mid_ind=round(length(sol.y(6,:))/2);
    d_T0_bar_mid{Ra_ind}=sol.y(6,mid_ind);
    w_hat_rms_mid{Ra_ind}=data_w_hat_rms{Ra_ind}.y(mid_ind);
    T_hat_rms_mid{Ra_ind}=data_T_hat_rms{Ra_ind}.y(mid_ind);
    u_hat_rms_mid{Ra_ind}=data_u_hat_rms{Ra_ind}.y(mid_ind);
    Nu{Ra_ind}=-sol.y(6,1);
end
%%compare u rms
plot_config.label_list={1,'$z$','$w_{rms}$'};
plot_config.name=['C:\Figure\DDC_LST\darcy_hewitt_w_rms',darcy_hewitt,'.png'];
plot_line(data_w_hat_rms,plot_config);

%%compare w rms
plot_config.label_list={1,'$z$','$u_{rms}$'};
plot_config.name=['C:\Figure\DDC_LST\darcy_hewitt_u_rms',darcy_hewitt,'.png'];
plot_line(data_u_hat_rms,plot_config);

%%compare T rms
plot_config.label_list={1,'$z$','$T_{rms}$'};
plot_config.name=['C:\Figure\DDC_LST\darcy_hewitt_T_rms',darcy_hewitt,'.png'];
plot_line(data_T_hat_rms,plot_config);

%%compare averaged temperature gradient
plot_config.label_list={1,'$\langle\bar{T}\rangle$','$z$'};
plot_config.fontsize_legend=24;
plot_config.print_size=[1,600,1000];
plot_config.name=['C:\Figure\DDC_LST\darcy_hewitt_T0_bar',darcy_hewitt,'.png'];
plot_line(data_T0_bar,plot_config);
plot_config.xlim_list=[1,0.45,0.55];
plot_config.name=['C:\Figure\DDC_LST\darcy_hewitt_T0_bar_zoom_in',darcy_hewitt,'.png'];
plot_line(data_T0_bar,plot_config);

function dz_dtheta = bvpfcn_2D(z,theta)
global Ra kx;
dz_dtheta = zeros(size(theta));
dz_dtheta = [theta(2)
             kx^2*theta(1)-kx^2*theta(3)
            theta(4)
            Ra*theta(1)*theta(6)+kx^2*theta(3)
            theta(6)
            2*Ra*theta(1)*theta(4)+2*Ra*theta(3)*theta(2)
            ];
end

function dz_dtheta = bvpfcn_3D(z,theta)
global Ra kx;
dz_dtheta = zeros(size(theta));
dz_dtheta = [theta(2)
             2*kx^2*theta(1)-2*kx^2*theta(3)
            theta(4)
            Ra*theta(1)*theta(6)+2*kx^2*theta(3)
            theta(6)
            2*Ra*theta(1)*theta(4)+2*Ra*theta(3)*theta(2)
            ];
end


function res = bcfcn(ya,yb)
res = [ya(1)
       yb(1)
       ya(3)
       yb(3)
       ya(5)-1
       yb(5)];
end



function dz_dtheta = bvpfcn_2D_two_mode(z,theta)
global Ra kx kxpp;
dz_dtheta = zeros(size(theta));
dz_dtheta = [theta(2);
             kx^2*theta(1)-kx^2*theta(3);
            theta(4);
            Ra*theta(1)*theta(6)+kx^2*theta(3);
            theta(6);
            2*Ra*theta(1)*theta(4)+2*Ra*theta(3)*theta(2)+2*Ra*theta(7)*theta(10)+2*Ra*theta(9)*theta(8);
            theta(8);
            kxpp^2*theta(7)-kxpp^2*theta(9);
            theta(10);
            Ra*theta(7)*theta(6)+kxpp^2*theta(9)];
end

function dz_dtheta = bvpfcn_3D_two_mode(z,theta)
global Ra kx kxpp;
dz_dtheta = zeros(size(theta));
dz_dtheta = [theta(2);
             2*kx^2*theta(1)-2*kx^2*theta(3);
            theta(4);
            Ra*theta(1)*theta(6)+2*kx^2*theta(3);
            theta(6);
            2*Ra*theta(1)*theta(4)+2*Ra*theta(3)*theta(2)+2*Ra*theta(7)*theta(10)+2*Ra*theta(9)*theta(8);
            theta(8);
            2*kxpp^2*theta(7)-2*kxpp^2*theta(9);
            theta(10);
            Ra*theta(7)*theta(6)+2*kxpp^2*theta(9)];
end


function res = bcfcn_two_mode(ya,yb)
res = [ya(1)
       yb(1)
       ya(3)
       yb(3)
       ya(5)-1
       yb(5)
       ya(7)
       yb(7)
       ya(9)
       yb(9)];
end

function res = bcfcn_periodic(ya,yb)
res = [ya(1)-yb(1)
       ya(2)-yb(2)
       ya(3)-yb(3)
       ya(4)-yb(4)
       ya(5)-yb(5)+1
       ya(6)-yb(6)];
end


function g = guess(x_num)
% for z_ind=1:length(d_theta_bar_num)
%    theta_num(z_ind)=sum(d_theta_bar_num(1:z_ind)); 
% end
global Ra kx darcy_hewitt
% g = [sin(pi*x_num)
%     pi*cos(pi*x_num)
%     sin(pi*x_num)
%     pi*cos(pi*x_num)
%     1-x_num-sin(2*pi*x_num)
%     -1-2*pi*cos(2*pi*x_num)
%      ];
W0=0.1;
switch darcy_hewitt
    case '2D'
        g=[W0*sin(pi*x_num)
            W0*pi*cos(pi*x_num)
            W0*(1+pi^2/(kx^2))*sin(pi*x_num)
            W0*pi*(1+pi^2/(kx^2))*cos(pi*x_num)
            1-x_num
            -1];
    case '3D'
        g=[W0*sin(pi*x_num)
            W0*pi*cos(pi*x_num)
            W0*(1+pi^2/(2*kx^2))*sin(pi*x_num)
            W0*pi*(1+pi^2/(2*kx^2))*cos(pi*x_num)
            1-x_num
            -1];
end
end


function g = guess_two_mode(x_num)
% for z_ind=1:length(d_theta_bar_num)
%    theta_num(z_ind)=sum(d_theta_bar_num(1:z_ind)); 
% end
global Ra kx darcy_hewitt kxpp
W0=0.1; W1=0.05;
switch darcy_hewitt
    case '2D'
        g=[W0*sin(pi*x_num)
            W0*pi*cos(pi*x_num)
            W0*(1+pi^2/(kx^2))*sin(pi*x_num)
            W0*pi*(1+pi^2/(kx^2))*cos(pi*x_num)
            1-x_num
            -1
            W1*sin(pi*x_num)
            W1*pi*cos(pi*x_num)
            W1*(1+pi^2/(kxpp^2))*sin(pi*x_num)
            W1*pi*(1+pi^2/(kxpp^2))*cos(pi*x_num)];
    case '3D'
        g=[W0*sin(pi*x_num)
            W0*pi*cos(pi*x_num)
            W0*(1+pi^2/(2*kx^2))*sin(pi*x_num)
            W0*pi*(1+pi^2/(2*kx^2))*cos(pi*x_num)
            1-x_num
            -1
            W1*sin(3*pi*x_num)
            W1*3*pi*cos(3*pi*x_num)
            W1*(1+pi^2/(2*kxpp^2))*sin(3*pi*x_num)
            W1*3*pi*(1+pi^2/(2*kxpp^2))*cos(3*pi*x_num)];
end

end


