clear all;
close all;
% 
% slurm_num={'12516590',
%     '12516677',
%     '12516679',
%     '12519906',
%     '12520491',
%     '12632332',
%     '12632508',
%     '12632582',
%     '12635575',
%     '12639413'};

%All of these are for porous media
%These are results comparing only for the thermal convection
% slurm_num={'12639319',
%         '12639413',
%         '12640102',
%         '12640104',
%         '12640192',
%         '12640193'};

%These are doing results for double-diffusive convection... 
% slurm_num={'12643033',
%             '12643076',
%             '12643064',
%             '12643209',
%             '12643124',
%             '12643067',
%             '12643068'};
% These for salt-finger regime...
% slurm_num={'12644479',
%            '12646311',
%            '12646317',
%            '12646319',
%            '12646320',
%            '12646324'}; %           '12646314',

slurm_num={'12678026'}; 
%'12640113', %%This case is wierd... not right
    
for slurm_ind=1:length(slurm_num)
    h5_name=['C:\Data\dedalus\dedalus_',...
        slurm_num{slurm_ind},...
            '\analysis\analysis_s1.h5'];
    dedalus_post_my{slurm_ind}=dedalus_post(h5_name);
    dedalus_post_my{slurm_ind}.uvw_hewitt=0;
    dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.dedalus_post_bvp();
    dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.bvp_plot;
        
end


plot_config.name=['C:\Figure\DDC_LST\HB_porous_','d_T_0.png'];
plot_line(data,plot_config)
% 
% for slurm_ind=length(slurm_num)
%     h5_name=['C:\Data\dedalus\dedalus_',...
%             slurm_num{slurm_ind},...
%                 '\analysis\analysis_s1.h5'];
% %             '\'];
% %     files = dir(h5_name);
% %     for folder_ind=1:length(files)
% %         h5_name=
%             
%     %     h5_name=['C:\Data\dedalus\dedalus_',...
%     %             slurm_num{slurm_ind},...
%     %             '\data.h5'];
%         h5disp(h5_name);
%         obj.z_list=h5read(h5_name,'/scales/z/1.0');
%         obj.w_hat=h5read(h5_name,'/tasks/w_hat');
%         obj.p_hat=h5read(h5_name,'/tasks/p_hat');
%         obj.T_hat=h5read(h5_name,'/tasks/T_hat');
%         obj.d_T_hat=h5read(h5_name,'/tasks/d_T_hat');
%         obj.S_hat=h5read(h5_name,'/tasks/S_hat');
%         obj.d_S_hat=h5read(h5_name,'/tasks/d_S_hat');
%         try
%             obj.w_hat_2=h5read(h5_name,'/tasks/w_hat_2');
%             obj.p_hat_2=h5read(h5_name,'/tasks/p_hat_2');
%             obj.T_hat_2=h5read(h5_name,'/tasks/T_hat_2');
%             obj.d_T_hat_2=h5read(h5_name,'/tasks/d_T_hat_2');
%             obj.S_hat_2=h5read(h5_name,'/tasks/S_hat_2');
%             obj.d_S_hat_2=h5read(h5_name,'/tasks/d_S_hat_2');
%         catch
%             disp('No second harmonic')
%         end
%         obj.T_0=h5read(h5_name,'/tasks/T_0');
%         obj.d_T_0=h5read(h5_name,'/tasks/d_T_0');
%         obj.S_0=h5read(h5_name,'/tasks/S_0');
%         obj.d_S_0=h5read(h5_name,'/tasks/d_S_0');
% 
%         obj.dy_T_mean=-1;
%         obj.dy_S_mean=-1;
%         data{1}.x=obj.T_0+1+obj.dy_T_mean*obj.z_list;
%         data{1}.y=obj.z_list;
%         plot_config.fontsize=20;
%         plot_config.label_list={1,'$\bar{T}_0+1+\bar{\mathcal{T}}_z z$', '$z$'};
%         plot_config.print_size=[1,500,900];
%         plot_config.name=['C:\Figure\DDC_LST\HB_porous_','T_0.png'];
%         plot_line(data,plot_config);
% 
%         data{1}.x=obj.d_T_0+obj.dy_T_mean;
%         data{1}.y=obj.z_list;
%         plot_config.fontsize=20;
%         plot_config.label_list={1,'$\partial_z \bar{T}_0+\bar{\mathcal{T}}_z$', '$z$'};
%         plot_config.print_size=[1,500,900];
%         plot_config.name=['C:\Figure\DDC_LST\HB_porous_','d_T_0.png'];
%         plot_line(data,plot_config);
%         
%         data{1}.x=obj.S_0+1+obj.dy_S_mean*obj.z_list;
%         data{1}.y=obj.z_list;
%         plot_config.label_list={1,'$\bar{S}_0+1+\bar{\mathcal{S}}_z z$', '$z$'};
%         plot_config.print_size=[1,500,900];
%         plot_config.name=['C:\Figure\DDC_LST\HB_porous_','S_0.png'];
%         plot_line(data,plot_config);
% 
%         data{1}.x=obj.S_0+1+obj.dy_S_mean*obj.z_list;
%         data{1}.y=obj.z_list;
%         plot_config.label_list={1,'$\bar{S}_0+1+\bar{\mathcal{S}}_z z$', '$z$'};
%         plot_config.print_size=[1,500,900];
%         plot_config.ylim_list=[1,0,0.01];
%         plot_config.name=['C:\Figure\DDC_LST\HB_porous_','S_0_local.png'];
%         plot_line(data,plot_config);
%         plot_config.ylim_list=0;
% 
%         data{1}.x=obj.S_0+1+obj.dy_S_mean*obj.z_list;
%         data{1}.y=obj.z_list;
%         plot_config.label_list={1,'$\bar{S}_0+1+\bar{\mathcal{S}}_z z$', '$z$'};
%         plot_config.print_size=[1,500,900];
%         plot_config.ylim_list=[1,0.49,0.51];
%         plot_config.xlim_list=[1,0.49,0.51];
%         plot_config.name=['C:\Figure\DDC_LST\HB_porous_','S_0_local_core.png'];
%         plot_line(data,plot_config);
%         plot_config.ylim_list=0; plot_config.xlim_list=0;
% 
% 
%         data{1}.x=obj.T_hat;
%         data{1}.y=obj.z_list;
%         plot_config.label_list={1,'$\widehat{T}$', '$z$'};
%         plot_config.print_size=[1,500,900];
%         plot_config.name=['C:\Figure\DDC_LST\HB_porous_','T_hat.png'];
%         plot_line(data,plot_config);
% 
%         data{1}.x=obj.S_hat;
%         data{1}.y=obj.z_list;
%         plot_config.label_list={1,'$\widehat{S}$', '$z$'};
%         plot_config.print_size=[1,500,900];
%         plot_config.name=['C:\Figure\DDC_LST\HB_porous_','S_hat.png'];
%         plot_line(data,plot_config);
% 
% 
%         data{1}.x=obj.S_hat;
%         data{1}.y=obj.z_list;
%         plot_config.label_list={1,'$\widehat{S}$', '$z$'};
%         plot_config.print_size=[1,500,900];
%         plot_config.ylim_list=[1,0,0.01];
%         plot_config.name=['C:\Figure\DDC_LST\HB_porous_','S_hat_local.png'];
%         plot_line(data,plot_config);
%         plot_config.ylim_list=0;
% 
%         data{1}.x=obj.w_hat;
%         data{1}.y=obj.z_list;
%         plot_config.label_list={1,'$\widehat{w}$', '$z$'};
%         plot_config.print_size=[1,500,900];
%         plot_config.name=['C:\Figure\DDC_LST\HB_porous_','w_hat.png'];
%         plot_line(data,plot_config);
% %     end
% end


