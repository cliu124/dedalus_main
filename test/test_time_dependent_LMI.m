clear all;
close all;
clc;


% A=[-1-9*cos(6*t)^2+6*sin(12*t), 12*cos(6*t)^2+9/2*sin(12*t);
%     -12*sin(6*t)^2+9/2*sin(12*t), -1-9*sin(6*t)^2-6*sin(12*t)];
t_list=linspace(0,2*pi,10);
g=9.8;
% l=sin(t_list);
for t_ind=1:length(t_list)
    t=t_list(t_ind);
    A_time=[-1-9*cos(6*t)^2+6*sin(12*t), 12*cos(6*t)^2+9/2*sin(12*t);
     -12*sin(6*t)^2+9/2*sin(12*t), -1-9*sin(6*t)^2-6*sin(12*t)];
%     l=10-sin(t);
%     A_time=[0,1;
%             -g/l,-1];
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

