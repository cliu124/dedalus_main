clear all;
close all;
clc;

syms kx ky kz Pr Ra_Tq nabla_2 nabla_p_2

A=[nabla_p_2/nabla_2*Pr*Ra_Tq, Pr*nabla_2;
    nabla_2,1];
A_eig=simplify(eig(A));


