clear all;
close all;
clc;


% A=[-1-9*cos(6*t)^2+6*sin(12*t), 12*cos(6*t)^2+9/2*sin(12*t);
%     -12*sin(6*t)^2+9/2*sin(12*t), -1-9*sin(6*t)^2-6*sin(12*t)];
t_list=linspace(0,2*pi,10);
g=9.8;
% l=sin(t_list);
T_list=[20:10:100];
Amp_list=[0.13:-0.01:0.07]/2;

for Amp_ind=1:length(Amp_list)
    for T_ind=1:length(T_list)
        for t_ind=1:length(t_list)
            t=t_list(t_ind);
        %     A_time=[-1-9*cos(6*t)^2+6*sin(12*t), 12*cos(6*t)^2+9/2*sin(12*t);
        %      -12*sin(6*t)^2+9/2*sin(12*t), -1-9*sin(6*t)^2-6*sin(12*t)];
            T=T_list(T_ind); 
            Amp=Amp_list(Amp_ind);
            l=1-Amp*sin(2*pi/T*t); %T=2, omega=2pi/T
            c=0;
            A_time=[0,1;
                    -g/l,-c];
            A11(t_ind,1)=A_time(1,1);
            A12(t_ind,1)=A_time(1,2);
            A{t_ind}=A_time;
        end


        sdpvar t

        P = sdpvar(2,2);
        Constraints = [P>=eye(2)];
        for A_ind=1:length(A)
            Constraints=[Constraints,A{A_ind}'*P+P*A{A_ind}+2*t*P<=0]
        end

        Objective = -t;
        diagnostics = bisection(Constraints, Objective,sdpsettings('solver','bisection','bisection.solver','sedumi','debug',1))

        alpha_list(Amp_ind,T_ind)=-value(t);
    end
end

save('alpha_list.mat','alpha_list');
data{1}.x=Amp_list;
data{1}.y=alpha_list(:,1);
plot_config.label_list={1,'A','$\lambda$'};
plot_config.name='time_varying_lambda_A_1.png';
plot_line(data,plot_config);

plot_config.legend_list={1,'$T=20s$'};
plot_config.name='time_varying_lambda_A_1_legend.png';
plot_line(data,plot_config);

data{1}.x=Amp_list;
data{1}.y=alpha_list(:,1);
data{2}.x=Amp_list;
data{2}.y=alpha_list(:,2);
plot_config.label_list={1,'A','$\lambda$'};
plot_config.name='time_varying_lambda_A_2.png';
plot_line(data,plot_config);

plot_config.legend_list={1,'$T=20s$','$T=40s$'};
plot_config.name='time_varying_lambda_A_2_legend.png';
plot_line(data,plot_config);


plot_config.legend_list={1};
for T_ind=1:length(T_list)
   data{T_ind}.x=Amp_list;
   data{T_ind}.y=alpha_list(:,T_ind);
   plot_config.legend_list{T_ind+1}=['T=',num2str(T_list(T_ind))];
end
plot_config.label_list={1,'A','$\lambda$'};
plot_config.Markerindex=3;
plot_config.user_color_style_marker_list=...
    {'k-','b-','r-','k--','b--','r--','k-.','b-.','r-.'};
plot_config.fontsize_legend=16;
plot_config.xlim_list=[0,0,0.07];
plot_config.name='time_varying_lambda_A_all_legend.png';
plot_line(data,plot_config);

clear data plot_config
plot_config.legend_list={1};
for Amp_ind=1:length(Amp_list)
   data{Amp_ind}.x=T_list;
   data{Amp_ind}.y=alpha_list(Amp_ind,:);
   plot_config.legend_list{Amp_ind+1}=['A=',num2str(Amp_list(Amp_ind))];
end
plot_config.label_list={1,'T','$\lambda$'};
plot_config.Markerindex=3;
plot_config.user_color_style_marker_list=...
    {'k-','b-','r-','k--','b--','r--','k-.','b-.','r-.'};
plot_config.fontsize_legend=16;
plot_config.name='time_varying_lambda_T_all_legend.png';
plot_line(data,plot_config);
