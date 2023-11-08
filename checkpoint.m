clear all;
close all;

slurm_num='4057882';
read_folder=['D:\Data\dedalus\dedalus_',slurm_num,'\analysis'];
write_folder=['D:\Data\dedalus\dedalus_',slurm_num,'\checkpoint'];

T=h5read_complex([read_folder,'\analysis_s1.h5'],'/tasks/T'); 
T_new=T(:,:,end);
h5write([write_folder,'\checkpoint_s1.h5'],'/tasks/T',T_new);

w=h5read_complex([read_folder,'\analysis_s1.h5'],'/tasks/w'); 
w_new=w(:,:,end);
h5write([write_folder,'\checkpoint_s1.h5'],'/tasks/w',w_new);

u=h5read_complex([read_folder,'\analysis_s1.h5'],'/tasks/u'); 
u_new=u(:,:,end);
h5write([write_folder,'\checkpoint_s1.h5'],'/tasks/u',u_new);

error('done');
p=h5read_complex([read_folder,'\analysis_s1.h5'],'/tasks/p'); 
p_new=p(:,:,end);
h5write([write_folder,'\checkpoint_s1.h5'],'/tasks/p',p_new);

