clear all;
close all;

global C R0 R2 theta_10_fun d_theta_10_fun d_theta_bar_fun kx;
flow='darcy_hewitt';
switch flow
    case 'benard'
        C=3;
        R2_list=[30,35,40,50,60,65,100,150];
        R0_list=[1,1,1,1,1,1];
        data_A{2}.y=1/24*R2_list.^2-R2_list;
        plot_config.legend_list={1,'$R_2=30$','$R_2=35$','$R_2=40$','$R_2=50$',...
            '$R_2=60$','$R_2=65$','$R_2=100$','$R_2=150$'};
        title_high_R2='$\theta_{10}=\frac{R_2}{2\sqrt{3}}[tanh\frac{1}{12}R_2z+tanh\frac{1}{12}R_2(1-z)-1]$';
        legend_A_2='$A_2=\frac{1}{24}R_2^2-R_2$';
    case 'benard_rotation'
%         kx=1;
%         R2_list=[10.89,20,40,60,100];
%         kx=0.5;
%         R2_list=[40,60,80,100,200,300];
        kx=2;
        R2_list=[20,40,60,80];
        data_A{2}.y=0*R2_list;
        plot_config.legend_list={1,'$R_0=10.89$','$R_0=20$','$R_0=40$','$R_0=60$',...
            '$R_0=100$'};
        title_high_R2='';
        legend_A_2='';
        %title_high_R2='$\theta_{10}=\frac{R_2}{2\sqrt{3}}[tanh\frac{1}{12}R_2z+tanh\frac{1}{12}R_2(1-z)-1]$';
        %legend_A_2='$A_2=\frac{1}{24}R_2^2-R_2$';
    case 'darcy'
        C=2;
        R2_list=[20,25,30,40,50,70];
        R0_list=[1,1,1,1,1,1];
        data_A{2}.y=1/16*R2_list.^2-R2_list;
        plot_config.legend_list={1,'$R_2=20$','$R_2=25$','$R_2=30$','$R_2=40$',...
            '$R_2=50$','$R_2=70$'};
        title_high_R2='$\theta_{10}=\frac{R_2}{2\sqrt{2}}[tanh\frac{1}{8}R_2z+tanh\frac{1}{8}R_2(1-z)-1]$';
        legend_A_2='$A_2=\frac{1}{16}R_2^2-R_2$';
    case 'darcy_hewitt'
        C=2;
        darcy_hewitt='2D';
        switch darcy_hewitt
            case '2D'
                Ra=10000;%[10000,20000,40000];
                kx_list=0.48*Ra.^0.4;
%                 R0_list=(1/0.48)^(1/0.4);
                R2_final=6000;
                R0_final=(Ra-R2_final)/kx_list^2;
                N=6;
%                 R2_list=R2_final*ones(1,6);
%                 R0_list=linspace(1,R0_final,N);
                R2_list=R2_final;
                R0_list=R0_final;
%                 R2_list=[10,10,10];
%                 R0_list=(Ra-R2_list)./kx_list.^2;
%                 R1_list=(Ra-kx_list.^2)./kx_list;
            case '3D'
                Ra=[4000,8000,16000];
                kx_list=0.17*Ra.^0.52;
%                 R2_list=[0,0,0];%Ra-R0_list*kx_list.^2;
%                 R0_list=(Ra-R2_list)./kx_list.^2;%(1/0.17)^(1/0.52);
%                 R1_list=(Ra-kx_list.^2)./kx_list;
%                 
                Ra=16000;
                kx_list=0.17*Ra^0.52;
                R2_list=Ra-kx_list^2;
                R0_list=1;%(Ra-R2_list)/kx_list^2;
                R1_list=(Ra-kx_list.^2)/kx_list;
        end
        data_A{2}.y=1/16*R2_list.^2-R2_list;

%         R2_list=Ra-kx_list.^2;
        %epsilon_list=kx_list.^(-2);
        %data_A{2}.y=1/16*R2_list.^2-R2_list;
        %plot_config.legend_list={0,'$R_2=20$','$R_2=25$','$R_2=30$','$R_2=40$',...
        %    '$R_2=50$','$R_2=70$'};
        title_high_R2='';
        legend_A_2='';
    
end
data_A{2}.x=R2_list;


% C=3; %This is 2 for porous media, and 3 for Benard problem as Blennerhassett & Bassom (1994)
% R2=150; %This is the parameter as eigenvalue...
for R2_ind=1:length(R2_list)
    R2=R2_list(R2_ind);
    R0=R0_list(R2_ind);
    syms x;
    switch flow
        case 'benard'
            %setup initial guess. This is the high R_2 approximation in
            %Blennerhassett & Bassom (1994)
            A2=1/24*R2^2-R2;
            theta_10=R2/2/sqrt(3)*(-1+tanh(1/12*R2*x)+tanh(1/12*R2*(1-x)));
            d_theta_bar_fun=matlabFunction(1/2*theta_10^2-A2);
        case {'darcy','darcy_hewitt'}
            %setup initial guess. This is the high R_2 approximation in
            %Lewis, Rees, Bassom (1997)
            A2=1/16*R2^2-R2;
            theta_10=R2/2/sqrt(2)*(-1+tanh(1/8*R2*x)+tanh(1/8*R2*(1-x)));
            d_theta_bar_fun=matlabFunction(1/2*theta_10^2-A2);

        case 'benard_rotation'
            %This seems to be a good initial guess
            theta_10=kx*R2/3*sin(pi*x);
            d_theta_bar=1-(-R2/3/(-R2*kx^2*R2/3)-kx^6/(-R2*kx^2));
            d_theta_bar_fun=matlabFunction(d_theta_bar,'Vars',x);
    end

    theta_10_fun=matlabFunction(theta_10);
    d_theta_10_fun=matlabFunction(diff(theta_10,x));

    xmesh=linspace(0,1,200);
    solinit=bvpinit(xmesh,@guess);
    switch flow
        case {'benard','darcy','darcy_hewitt'}
            sol=bvp4c(@bvpfcn,@bcfcn,solinit);
%             sol=bvp4c(@bvpfcn,@bc_periodic_fcn,solinit);
        case 'benard_rotation'
            sol=bvp4c(@bvpfcn_rotation,@bcfcn,solinit);
    end
    A_num(R2_ind,1)=1/2*sum(sol.y(3,1:end-1).^2.*diff(sol.x));
    

    %syms x;
    %theta_10=R2/2/sqrt(3)*(-1+tanh(1/12*R2*x)+tanh(1/12*R2*(1-x)));
    %res_theta_10=3*diff(theta_10,x,2)+R2*theta_10-1/2*theta_10*(theta_10^2-2*(1/24*R2^2-R2));
    %dz=diff(xmesh); dz=dz(1);
    %A_ana=1/2*sum(double(subs(theta_10,x,xmesh)).^2)*dz;

    %res=A-(1/24*R2^2-R2);
    
    data_theta_bar{R2_ind}.x=sol.x;
    data_theta_bar{R2_ind}.y=sol.y(1,:);
    
    data_d_theta_bar{R2_ind}.x=sol.x;
    data_d_theta_bar{R2_ind}.y=sol.y(2,:);
    
    data_theta_10{R2_ind}.x=sol.x;
    data_theta_10{R2_ind}.y=sol.y(3,:);
    
    data_theta_10_high_R2{R2_ind}.x=sol.x;
    switch flow
        case {'benard','benard_rotation'}
            data_theta_10_high_R2{R2_ind}.y=...
                R2/2/sqrt(3)*(-1+tanh(1/12*R2*sol.x)+tanh(1/12*R2*(1-sol.x)));
        case {'darcy','darcy_hewitt'}
            data_theta_10_high_R2{R2_ind}.y=...
                R2/2/sqrt(2)*(-1+tanh(1/8*R2*sol.x)+tanh(1/8*R2*(1-sol.x)));
    end
    data_A{1}.x(R2_ind)=R2;
    data_A{1}.y(R2_ind)=A_num(R2_ind,1);
end
d_theta_bar=data_d_theta_bar{end}.y(1,round(length(data_d_theta_bar{end}.y)/2));
d_theta_bar_rescale=d_theta_bar*kx_list^(-2);
theta_bar_full=1-data_theta_bar{end}.x+data_theta_bar{end}.y*kx_list^(-2);
plot(data_theta_bar{end}.x,theta_bar_full);


switch flow
    case {'benard','darcy'}
        plot_config.label_list={1,'$z$','$\bar{\theta}_2$'};
        plot_config.Markerindex=3;
        plot_config.ytick_list=[1,-50,-40,-30,-20,-10,0,10,20,30,40,50];
        plot_config.xtick_list=[1,0,0.2,0.4,0.6,0.8,1];
        plot_config.user_color_style_marker_list={'k-','k--','b-','b--','r-','r--','m-','m--'};
        plot_config.fontsize_legend=20;
        plot_config.name=['C:\Figure\DDC_LST\',flow,'_theta_bar.png'];
        plot_line(data_theta_bar,plot_config);

        %%plot of theta_10
        plot_config.label_list={1,'$z$','$\theta_{10}$'};
        plot_config.name=['C:\Figure\DDC_LST\',flow,'_theta_10.png'];
        plot_line(data_theta_10,plot_config);

        plot_config.label_list={1,'$z$','$\theta_{10}$'};
        plot_config.title_list={1,title_high_R2};
        plot_config.name=['C:\Figure\DDC_LST\',flow,'_theta_10_high_R2.png'];
        plot_line(data_theta_10_high_R2,plot_config);

        %plot of A over the analytical expression...

        plot_config.title_list={0};
        plot_config.name=['C:\Figure\DDC_LST\',flow,'_A.png'];
        plot_config.label_list={1,'$R_2$','$\bar{A}_2$'};
        plot_config.user_color_style_marker_list={'k*','b-'};
        plot_config.ytick_list=[0,0,20,40,60,80,100,120,150];
        plot_config.xtick_list=[1,30,50,70,90,110,130,150];
        plot_config.legend_list={1,'Numerical',legend_A_2};
        plot_config.fontsize_legend=30;
        plot_line(data_A,plot_config);
    case 'darcy_hewitt'
        
        plot_config.label_list={1,'$z$','$\bar{\theta}_2$'};
        plot_config.Markerindex=3;
        plot_config.user_color_style_marker_list={'k-','k--','b-','b--','r-','r--','m-','m--'};
        plot_config.fontsize_legend=20;
        plot_config.name=['C:\Figure\DDC_LST\',flow,'_theta_bar.png'];
        plot_line(data_theta_bar,plot_config);

        %%plot of total temperature
        
        data_theta_bar_total=data_theta_bar;
        data_theta_bar_total{1}.y=1-data_theta_bar{1}.x+data_theta_bar{1}.y*kx_list^(-2);
        plot_config.label_list={1,'$z$','$1-z+k_x^{-2}\bar{\theta}_2$'};
        plot_config.Markerindex=3;
        plot_config.user_color_style_marker_list={'k-','k--','b-','b--','r-','r--','m-','m--'};
        plot_config.fontsize_legend=20;
        plot_config.name=['C:\Figure\DDC_LST\',flow,'_theta_bar_total.png'];
        plot_line(data_theta_bar_total,plot_config);

        
        %%plot of theta_10
        plot_config.label_list={1,'$z$','$\theta_{10}$'};
        plot_config.name=['C:\Figure\DDC_LST\',flow,'_theta_10.png'];
        plot_line(data_theta_10,plot_config);
        
        
    case 'benard_rotation'
        plot_config.label_list={1,'$z$','$\bar{\theta}_0$'};
        plot_config.Markerindex=3;
        plot_config.ytick_list=[1,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1];
        plot_config.xtick_list=[1,0,0.2,0.4,0.6,0.8,1];
        plot_config.user_color_style_marker_list={'k-','k--','b-','b--','r-','r--','m-','m--'};
        plot_config.fontsize_legend=20;
        plot_config.name=['C:\Figure\DDC_LST\',flow,'_theta_bar.png'];
        for R2_ind=1:length(R2_list)
            data_theta_bar{R2_ind}.y=1-data_theta_bar{R2_ind}.x+data_theta_bar{R2_ind}.y;
            data_theta_10{R2_ind}.y=2*data_theta_10{R2_ind}.y*kx;
        end
        plot_line(data_theta_bar,plot_config);

        %%plot of theta_10      
        
        plot_config.ytick_list=[0,0,20,40,60,80,100,120,150];
        plot_config.label_list={1,'$z$','$W_{10}$'};
        plot_config.name=['C:\Figure\DDC_LST\',flow,'_W_10.png'];
        plot_line(data_theta_10,plot_config);

end

%figure(1)
%plot(sol.x,sol.y(1,:));
%figure(2)
%plot(sol.x,sol.y(3,:)); hold on;
%plot(sol.x,R2/2/sqrt(3)*(-1+tanh(1/12*R2*sol.x)+tanh(1/12*R2*(1-sol.x))),'r--');
% theta_10_sym=R2_sym/2/sqrt(3)*(-1+tanh(1/12*R2_sym*x)+tanh(1/12*R2_sym*(1-x)));
% res_theta_10_sym=3*diff(theta_10,x,2)+R2_sym*theta_10-1/2*theta_10*(theta_10^2-2*(1/24*R2_sym^2-R2_sym));

% figure(3);
% plot(xmesh,double(subs(res_theta_10,x,xmesh)));
%A_sym=int(1/2*theta_10^2,x,0,1);

function dz_dtheta = bvpfcn(z,theta)
global C R2 R0;
dz_dtheta = zeros(size(theta));
dz_dtheta = [theta(2)
            theta(3)*theta(4)
            theta(4)
            1/C*(-R2*theta(3)+R0*theta(3)*theta(2) )];
end

function dz_dtheta = bvpfcn_rotation(z,theta)
global R2 kx;
dz_dtheta = zeros(size(theta));
dz_dtheta = [theta(2)
    -4*theta(3)*theta(4)*theta(2)/(1+2*theta(3)^2)+4*theta(3)*theta(4)/(1+2*theta(3)^2)
    theta(4)
    kx^6*theta(3)-R2*kx^2*(1-theta(2))*theta(3)];

end

function res = bcfcn(ya,yb)
res = [ya(1)
       yb(1)
       ya(3)
       yb(3)];
end

function res = bc_periodic_fcn(ya,yb)
res = [ya(1)-yb(1)
       ya(2)-yb(2)
       ya(3)-yb(3)
       ya(4)-yb(4)];
end

function g = guess(x_num)
global R2 theta_10_fun d_theta_10_fun d_theta_bar_fun;;
% for z_ind=1:length(d_theta_bar_num)
%    theta_num(z_ind)=sum(d_theta_bar_num(1:z_ind)); 
% end

g = [-R2/2*sin(2*pi*x_num)
     d_theta_bar_fun(x_num)
     theta_10_fun(x_num)
     d_theta_10_fun(x_num)
     ];
end
