clear all;
close all;

slurm_num={'20230307231811',...
    '20230307232436',...
    '20230307232924'};
slurm_num={'20230320172744'};
% slurm_num={'20230308000459'};
% slurm_num={'20230308101259'};
slurm_num={'20230320172744'};
slurm_num={'20230314185834',
'20230314190117',
'20230314190403',
'20230314190655',
'20230314190928',
'20230314191201',
'20230314191934',
'20230314193155',
'20230314200855',
'20230314201731',
'20230314211543',
'20230314213422',
'20230314220008',
'20230314221032',
'20230314224402',
'20230314224941',
'20230314225825',
'20230314231644',
'20230314232329',
'20230314233351',
'20230315000200'
};
slurm_num={'20230326222249'};
flag.print=1;
flag.visible=0;
flag.video=0;
flag.no_ylabel=0;
flag.post_plot=1;

% for slurm_ind=1:length(slurm_num)
%     h5_name=['D:\Data\dedalus\dedalus_',...
%         slurm_num{slurm_ind},...
%         '\analysis\analysis_s1.h5'];
% 
%      set(0,'DefaultFigureVisible','on')
%      dedalus_post_my{slurm_ind}=dedalus_post(h5_name,flag);
% %      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.dedalus_post_ivp();
%      dedalus_post_my{slurm_ind}.print=1; dedalus_post_my{slurm_ind}.visible=1;
%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.spectrum_t('U_0',[0.25],[],[2]);
%      
%      %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('u');
%      data{1}.y(slurm_ind)=dedalus_post_my{slurm_ind}.freq_sort(1);
%      data{1}.x(slurm_ind)=dedalus_post_my{slurm_ind}.Ra_T;
%      
% end
% data{1}.y=data{1}.y*2*pi;
% plot_config.label_list={1,'$Ra_{T,q}$','$\omega$'};
% plot_config.name='RBC_Ra_global_SM_Ra_Tq_omega.png';
% plot_line(data,plot_config);
% error('1');

for slurm_ind=1:length(slurm_num)
    h5_name=['D:\Data\dedalus\dedalus_',...
        slurm_num{slurm_ind},...
        '\analysis\analysis_s1.h5'];

     set(0,'DefaultFigureVisible','on')
     dedalus_post_my{slurm_ind}=dedalus_post(h5_name,flag);
     
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('dy_T_mean_q');     
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.get_Nu('T',[]);
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.total_xt_ave('T',[],[]);
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('u');
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('T');
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('w');
     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.phase_diagram('u','T',[],'max_z2')

end