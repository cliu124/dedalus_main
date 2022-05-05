clear all;
close all;
clc;

slurm_num={'13625643'};
flag.print=1;
flag.video=0;
flag.visible=0;
flag.no_ylabel=0;
for slurm_ind=1:length(slurm_num)%:length(slurm_num)-1%[find(strcmp(slurm_num,'12247549'))]%slurm_ind=length(slurm_num)-2:length(slurm_num)-1
    
    %%change the path into D... just store data in the external disk...
    h5_name=['.\dedalus_',...
        slurm_num{slurm_ind},...
        '\analysis\analysis_s1.h5'];

     set(0,'DefaultFigureVisible','on')
     dedalus_post_my{slurm_ind}=dedalus_post(h5_name,flag);
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.dedalus_post_ivp();
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.get_Nu('S',170);
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.E_time('S',0);

     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('S');
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('u');
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.rms_xt('u');
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.total_xt_ave('S');
     
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.z_slice('S_tot',0.5);

     dedalus_post_my{slurm_ind}.print=0; dedalus_post_my{slurm_ind}.visible=0;
     dedalus_post_my{slurm_ind}.video=0;
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.snapshot('S_tot',1);
    
end
