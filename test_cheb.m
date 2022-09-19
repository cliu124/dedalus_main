clear all;
close all;
clc;

N=30;
M=2;
[x, DM] = chebdif(N, M);
D1=DM(:,:,1);
D2=DM(:,:,2);



