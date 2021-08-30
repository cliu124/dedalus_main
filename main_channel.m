clear all;
close all;
clc;

folder_name='C:\Data\dedalus\channel\';
file_name='channel_s1';
h5_name=[folder_name,file_name,'.h5'];
h5disp(h5_name);
x=h5read(h5_name,'/scales/x/1');
% z=h5read(h5_name,'/scales/z/1');
z=h5read(h5_name,'/scales/y/1');
t=h5read(h5_name,'/scales/sim_time');

Lx=max(x)-min(x);
Lz=max(z)-min(z);

% u=h5read(h5_name,'/tasks/u midplane');
% w=h5read(h5_name,'/tasks/w midplane');
u=h5read(h5_name,'/tasks/u');
w=h5read(h5_name,'/tasks/w');
% T=h5read(h5_name,'/tasks/T');
% S=h5read(h5_name,'/tasks/S');
