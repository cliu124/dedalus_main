clear all;
close all;

slurm_num={'20230307130932'};
flag.print=1;
flag.visible=0;
flag.video=0;
flag.no_ylabel=0;
flag.post_plot=1;
for slurm_ind=1:length(slurm_num)
    h5_name=['D:\Data\dedalus\dedalus_',...
        slurm_num{slurm_ind},...
        '\analysis\analysis_s1.h5'];

     set(0,'DefaultFigureVisible','on')
     dedalus_post_my{slurm_ind}=dedalus_post(h5_name,flag);
     
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('dy_T_mean_q');     
     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.get_Nu('T',[]);
     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.total_xt_ave('T',[],[]);
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('u');
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('T');
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('w');
     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.phase_diagram('u','T',[],'max_z2')

end