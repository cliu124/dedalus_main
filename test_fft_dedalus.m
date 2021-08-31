clear all;
close all;
clc;

%%12073090: 32*32 size, elevator mode, without shear
%%12073116: 96*32 size, elevator mode, without shear
%%12075689:  
%%12075690:

%%These are local folder
% folder_name='C:\Data\dedalus\IFSC_2D_without_shear\';
%folder_name='C:\Data\dedalus\IFSC_2D_with_shear\';
%folder_name='C:\Data\dedalus\dedalus_12073090\IFSC_2D_without_shear\';
folder_name='C:\Data\dedalus\dedalus_12073116\IFSC_2D_without_shear\';
% folder_name='C:\Data\dedalus\dedalus_12075689\IFSC_2D_with_shear\';

% folder_name='/rc_scratch/chli3324/dedalus_12073090/';
% file_name='IFSC_2D_without_shear_s1_random';
% file_name='IFSC_2D_without_shear_s1_small_domain';

% file_name='IFSC_2D_without_shear_s1_elevator_short';
% file_name='IFSC_2D_without_shear_s1_elevator_long';
% file_name='IFSC_2D_with_shear_s1';
file_name='IFSC_2D_without_shear_s1';
h5_name=[folder_name,file_name,'.h5'];
h5disp(h5_name);
x=h5read(h5_name,'/scales/x/1.0');
Nx=length(x);
z=h5read(h5_name,'/scales/z/1.0');
Nz=length(z);
t=h5read(h5_name,'/scales/sim_time');
kx_list=h5read(h5_name,'/scales/kx');
kz_list=h5read(h5_name,'/scales/kz');
Lx=max(x)-min(x)+x(2);
Lz=max(z)-min(z)+z(2);

% u=h5read(h5_name,'/tasks/u');
% w=h5read(h5_name,'/tasks/w');
% T=h5read(h5_name,'/tasks/T');
S=h5read(h5_name,'/tasks/S');
S_coeff=h5read(h5_name,'/tasks/S_coeff');
S_coeff=S_coeff.r+1i*S_coeff.i;
t_fft_ind=50
S_fft2=fft2(S(:,:,t_fft_ind))/Nx/Nz;

kx_list_fft=(0:(2*pi/Lx):(2*pi/Lx)*(Nx/2-1))';
kz_list_fft=(0:(2*pi/Lz):(2*pi/Lz)*(Nz/2-1))';
error_kx=kx_list_fft-kx_list;
error_kz=kz_list_fft-kz_list(1:Nz/2);
%%Then, these two are equivalent... 
mesh(kx_list_fft,kz_list_fft,abs(S_fft2(1:length(kz_list_fft),1:length(kx_list_fft))));
mesh(kx_list,kz_list(1:Nz/2),abs(S_coeff(1:Nz/2,:,t_fft_ind)));
error_S_fft=norm(S_fft2(1:length(kz_list_fft),1:length(kx_list_fft))-S_coeff(1:Nz/2,:,t_fft_ind));
