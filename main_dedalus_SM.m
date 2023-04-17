clear all;
close all;

slurm_num={'20230307231811',...
    '20230307232436',...
    '20230307232924'};
slurm_num={'20230320172744'};
% slurm_num={'20230308000459'};
% slurm_num={'20230308101259'};
slurm_num={'20230320172744'};
slurm_num={
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
'20230315000200',
'20230403232615',
'20230404164722',
'20230404181136',
'20230404181720',
'20230404182254',
'20230404183245'
};
% slurm_num=slurm_num(end-6:end-7);
% slurm_num={'15249573'};
% slurm_num={'15295279'}
% slurm_num=slurm_num(end);
% slurm_num={'20230326211527'};
% slurm_num={'15241562'};
% slurm_num={'15243701',
%             '15243704',
%             '15243718',
%             '15243757',
%             '15243771',
%             '15243776',
%             '15243778'};
% slurm_num={'15243781',
%         '15243790',
%         '15243799',
%         '15243818',
%         '15243854',
%         '15243857',
%         '15243865'};
%     
% slurm_num={'15245854',
% '15247183',
% '15247191',
% '15247196',
% '15247200',
% '15247257',
% '15247260',
% '15247261'
% };
% slurm_num={'15247257'};
% slurm_num={'15245786',
% '15247389',
% '15247390',
% '15247396',
% '15247401',
% '15247411',
% '15247418',
% '15247419'};
% slurm_num=slurm_num(5:end);
% slurm_num={'20230404181720'};
% slurm_num={'20230410120319'};
slurm_num={'20230416120557'};
flag.print=1;
flag.visible=1;
flag.video=0;
flag.no_ylabel=0;
flag.post_plot=1;

for slurm_ind=1:length(slurm_num)
    h5_name=['D:\Data\dedalus\dedalus_',...
        slurm_num{slurm_ind},...
        '\analysis\analysis_s1.h5'];

     set(0,'DefaultFigureVisible','on')
     dedalus_post_my{slurm_ind}=dedalus_post(h5_name,flag);
%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.dedalus_post_ivp();
%      dedalus_post_my{slurm_ind}.print=0;
%      dedalus_post_my{slurm_ind}.visible=1;
     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.spectrum_t('U_0',[0.25],[],[2]);
     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.elevator_growing('T');
     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('u',[]);
     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('T',[]);

     %data{1}.y(slurm_ind)=dedalus_post_my{slurm_ind}.freq_sort(1);
     %data{1}.y(slurm_ind)=2*pi/dedalus_post_my{slurm_ind}.period_t;
%      data{1}.x(slurm_ind)=dedalus_post_my{slurm_ind}.Ra_T;
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('dy_T_mean_q');     
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.get_Nu('T',[200]);
     %data_Nu{1}.y(slurm_ind)=dedalus_post_my{slurm_ind}.Nu;
     
end
error('1');
% data{2}.x=data{1}.x;
% data_Nu{1}.x=data{1}.x;
% data{2}.x=data{1}.x(11:end);
Ra_g=46761.08197624290000;
fit_ind=8:length(data{1}.x);
modelfun = @(b,x)(b(1)./(-log(b(2)-x)));
mdl = fitnlm(data{1}.x(fit_ind),data{1}.y(fit_ind),modelfun,[80,Ra_g]);

data{2}.x=data{1}.x(fit_ind);
data{2}.y=80./(-log(Ra_g-data{2}.x));
% factor=data{1}.y(end-1)/data{2}.y(end-1);
% data{2}.y=factor*data{2}.y;
plot_config.label_list={1,'$Ra_{T,q}$','$\omega$'};
plot_config.name='RBC_Ra_global_SM_Ra_Tq_omega.png';
plot_config.user_color_style_marker_list={'msquare','k--'};
plot_config.Markerindex=3;
plot_config.print=1;
plot_config.xlim_list=[1,Ra_g-0.01,Ra_g];
plot_line(data,plot_config);
% 
% data{1}.x=data{1}.x(fit_ind);
% data{1}.y=data{1}.y(fit_ind);
% plot_config.user_color_style_marker_list={'msquare','k--'};

plot_line(data,plot_config);

error('1');

for slurm_ind=1:length(slurm_num)
    h5_name=['D:\Data\dedalus\dedalus_',...
        slurm_num{slurm_ind},...
        '\analysis\analysis_s1.h5'];

     set(0,'DefaultFigureVisible','on')
     dedalus_post_my{slurm_ind}=dedalus_post(h5_name,flag);
     
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('dy_T_mean_q');     
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.get_Nu('T',[1000]);
%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.total_xt_ave('T',[],[]);
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('u');
%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('T');
%      dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.x_ave('w');
     %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.phase_diagram('u','T',[],'max_z2')
     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.spectrum_t('u',[0.25],[],[2]);

end