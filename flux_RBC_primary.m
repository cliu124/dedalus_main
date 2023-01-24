clear all;
close all;
clc;
syms Ra_Tq Pr K k_perp;
A=[-Pr*K^2, Pr*Ra_Tq*k_perp^2/K^2;
    1, -K^2];
eig_A=eig(A);


