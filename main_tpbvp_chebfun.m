clear all;
close all;
%% System of two nonlinear BVPs
% Asgeir Birkisson and Toby Driscoll, September 2010

%%
% (Chebfun example ode-nonlin/BVPSystem.m)
% [Tags: #nonlinearODE, #ODEsystem]

%% System of equations
% Here is a system of two coupled nonlinear ODEs on the
% interval $[-1,1]$, with boundary conditions.
%
% $$ u'' - \sin(v) = 0, $$
%
% $$ v'' + \cos(u) = 0, $$               
%
% $$ u(-1) = 1,  ~~  v'(-1) = 0, ~~ u'(1) = 0,  ~~  v(1) = 0. $$

%% Solution using multiple variables `u` and `v`
% One way you can solve a problem like this with Chebfun is to work with
% multiple variables, solving for two chebfuns $u$ and $v$. Here we do this,
% setting up the problem using anonymous functions that take two chebfuns as
% input and return the two chebfuns as output:
% addpath('../chebfun/')
N = chebop(0, 1);
x = chebfun('x',[0,1]);
% R2=30; R0=3;
C=2;
N.lbc = @(u,v)[ u; v];
N.rbc =  @(u,v)[ u; v];
Ra=[4000];%,8000,16000
kx_list=0.17*Ra.^0.52;
%R2_list=[0];%Ra-R0_list*kx_list.^2;
%R0_list=(Ra-R2_list)./kx_list.^2;%(1/0.17)^(1/0.52);
R2_list=150;
R0_list=1;
R2_con_list=linspace(pi^2*C,R2_list,10);%%The list for performing continuation solving...
R0_con_list=linspace(0,R0_list,10);
% R2 = [1 .5 .2 .1 .03 .01 .003];
N.init = [-1/(8*pi)*sin(2*pi*x); sin(pi*x)];
for R2_con_ind=1:length(R2_con_list)
    for R0_con_ind=1:length(R0_con_list)
          R0=R0_con_list(R0_con_ind);
          R2=R2_con_list(R2_con_ind);
          N.op = @(x,u,v)[ 2*diff(u,2)+R2*u-R0*u*diff(v) ; diff(v,2)-4*u*diff(u) ];
%           N.init = 10*[sin(pi*x); sin(pi*x)];
          [u, v, info] = solvebvp(N,[0; 0]);
          N.init = [u; v];
          plot(u); hold on;
          plot(v);
          1
    end
end
nrmduvec = info.normDelta;

%%
% We can now plot the solution compomenents u and v:
LW = 'linewidth'; FS = 'fontsize';
figure, subplot(1,2,1), plot(u, LW, 2)
hold on, plot(v,'--r', LW, 2), hold off
title('u and v vs. x', FS, 10), legend('u', 'v')
box on, grid on
xlabel('x', FS, 10), ylabel('u(x) and v(x)', FS, 10)
subplot(1,2,2), semilogy(nrmduvec, '-*', LW, 2)
title('Norm of update vs. iteration no.', FS, 10)
box on, grid on
xlabel('iteration no.', FS, 10), ylabel('norm of update', FS, 10)

error('1')
%% Solution using a single indexed variable `u`
% Another way to solve the same problem is to work with a single chebmatrix
% variable `u` that has two components, `u{1}` and `u{2}`.
%
% $$ (u_1)'' - \sin(u_2) = 0, $$
%
% $$ (u_2)'' + \cos(u_1) = 0, $$               
%
% $$ u_1(-1) = 1, ~~ (u_2)'(-1) = 0, ~~ (u_1)'(1) = 0, ~~ u_2(1) = 0. $$
%
N = chebop(-1, 1);
x = chebfun('x');
N.op = @(x,u) [ diff(u{1},2) - sin(u{2}); diff(u{2},2) + cos(u{1}) ];
N.lbc = @(u)[ u{1} - 1; diff(u{2}) ];
N.rbc =  @(u)[ u{2}; diff(u{1}) ];
N.init = [0*x; 0*x];

%%
% The solution process is the same, but now the solutions gets returned as the
% chebmatrix u:
u = N\[0; 0];
nrmduvec = info.normDelta;

%%
% The now use the curly braces notation of chebmatrices to retrieve the solution
% componennts.
clf
plot(u{1}, LW, 2), hold on
plot(u{2}, '--r', LW, 2), hold off
title('u_1(x) and u_2(x) vs. x', FS, 10), legend('u_1', 'u_2')
box on, grid on
xlabel('x', FS, 10), ylabel('u_1(x) and u_2(x)', FS, 10)