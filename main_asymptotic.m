clear all;
close all;

clc;

syms x z t 
% independent=[x,z,t];
syms x_star z_star t_star;
syms epsilon
assume(epsilon,'real');
assume(epsilon,'positive');
syms R_0;
Ra=1+R_0*epsilon;%+R_0*epsilon^(-1); %%set up the Rayleigh number 

%%
x_x_star=1; %%x/x_star, x_star means dimensional quantity or before change of variable.... 
z_z_star=1;
t_t_star=1;

N=3;%set up the order of asymtotic expansion
syms psi_(x,z) [N,1]
syms T_(x,z) [N,1]
syms T_0(x,z) T_10(x,z)%set up additional zero order.
syms psi_0(x,z)
T=T_0+T_10+epsilon*T_1+epsilon^2*T_2+epsilon^3*T_3;
psi=psi_0+epsilon*psi_1+epsilon^2*psi_2+epsilon^3*psi_3;
%T=T_0+epsilon^2*T_1+epsilon^4*T_2+epsilon^6*T_3;
%psi=epsilon*psi_0+epsilon^3*psi_1+epsilon^5*psi_2+epsilon^7*psi_3;

%Porous media
flow='benard';%{'darcy','benard'}%
switch flow
    case 'darcy'
        equation1=diff(psi,x,2)*x_x_star^2+epsilon*diff(psi,z,2)*z_z_star^2+diff(T,x)*x_x_star*Ra;
    case 'benard'
        equation1=diff(psi,x,4)*x_x_star^2+2*epsilon*diff(diff(psi,x,2),z,2)+epsilon^2*diff(psi,z,4)*z_z_star^2-diff(T,x)*x_x_star*Ra;
end
equation2=epsilon*(diff(psi,z)*z_z_star*diff(T,x)*x_x_star-diff(psi,x)*x_x_star*diff(T,z)*z_z_star)-diff(T,x,2)*x_x_star^2-epsilon*diff(T,z,2)*z_z_star^2+diff(psi,x)*x_x_star;

%Benard problem
%equation1=diff(psi,x,2)*x_x_star^2+diff(psi,z,2)*z_z_star^2+diff(T,x)*x_x_star;
%equation2=diff(psi,z)*z_z_star*diff(T,x)*x_x_star-diff(psi,x)*x_x_star*diff(T,z)*z_z_star+diff(psi,x)*x_x_star-diff(T,x,2)*x_x_star^2/Ra-diff(T,z,2)*z_z_star^2/Ra;

exponent1=3;
exponent2=2;

expand1=formula(coeffs(simplify(expand(equation1)*epsilon^exponent1,50),epsilon));
expand2=formula(coeffs(simplify(expand(equation2)*epsilon^exponent2,50),epsilon));
expand1=expand1(1,1:N+1);
expand2=expand2(1,1:N+1);

latex_str1=['\begin{subequations}',newline,'\begin{align}',newline,...
    'O(\epsilon^{-3}):&\;\;',latex(expand1(1)),'=0\\',newline,...
    'O(\epsilon^{-2}):&\;\;',latex(expand1(2)),'=0\\',newline,...
    'O(\epsilon^{-1}):&\;\;',latex(expand1(3)),'=0\\',newline...
    'O(\epsilon^{0}):&\;\;',latex(expand1(4)),'=0',newline,...
    '\end{align}',newline,'\end{subequations}'];
latex_str2=['\begin{subequations}',newline,'\begin{align}',newline,...
    'O(\epsilon^{-2}):&\;\;',latex(expand2(1)),'=0\\',newline,...
    'O(\epsilon^{-1}):&\;\;',latex(expand2(2)),'=0\\',newline,...
    'O(\epsilon^{-0}):&\;\;',latex(expand2(3)),'=0\\',newline,...
    'O(\epsilon^1):&\;\;',latex(expand2(4)),'=0',newline,...
    '\end{align}',newline,'\end{subequations}'];

latex_str1=strrep(latex_str1,'\left(x,z\right)','');
latex_str1=strrep(latex_str1,'\frac{\partial }{\partial x}','\partial_x');
latex_str1=strrep(latex_str1,'\frac{\partial }{\partial z}','\partial_z');
latex_str1=strrep(latex_str1,'\frac{\partial ^2}{\partial x^2}','\partial_x^2');
latex_str1=strrep(latex_str1,'\frac{\partial ^2}{\partial z^2}','\partial_z^2');

latex_str2=strrep(latex_str2,'\left(x,z\right)','');
latex_str2=strrep(latex_str2,'\frac{\partial }{\partial x}','\partial_x');
latex_str2=strrep(latex_str2,'\frac{\partial }{\partial z}','\partial_z');
latex_str2=strrep(latex_str2,'\frac{\partial ^2}{\partial x^2}','\partial_x^2');
latex_str2=strrep(latex_str2,'\frac{\partial ^2}{\partial z^2}','\partial_z^2');


% 
% expand1=expand1(1,1:N+1);
% expand2=expand2(1,1:N+1);


% variable_list={'T','psi'};
% for variable_ind=1:length(variable_list)
%     for n=0:N
%         variable.([variable_list{variable_ind},num2str(n)])=syms([variable_list{variable_ind})
% 
%     end
% end




% T=syms('T',[N,1]);
