clear all;
close all;
clc;

sdpvar t
A = [-1 2;-3 -4];
P = sdpvar(2,2);
Constraints = [P>=eye(2), A'*P+P*A <= -2*t*P];
Objective = -t;
diagnostics = bisection(Constraints, Objective,sdpsettings('solver','bisection','bisection.solver','sedumi'))



%,sdpsettings('solver','mosek')

