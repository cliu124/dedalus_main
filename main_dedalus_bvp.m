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

group_name='HB_porous_thermal_BC';
%All of these are for porous media
switch group_name
    case 'HB_porous_thermal_BC'
        %These are results comparing only for the thermal convection
        slurm_num={'12639319',
                '12640193',
                '12640192',
                '12639413',
                '12640102'
                }; %'12640104',
    case 'HB_porous_diffusive_BC'
        %These are doing results for double-diffusive convection... 
        slurm_num={'12643033',
                    '12643076',
                    '12643064',
                    '12643209',
                    '12643124',
                    '12643067',
                    '12643068'};
    case 'HB_porous_finger_BC'
        % These for salt-finger regime...
        slurm_num={'12644479',
                   '12646311',
                   '12704259',
                   '12646317',
                   '12646319',
                   '12646320',
                   '12646324'}; %           '12646314',
    case 'hewitt_2D'
        slurm_num={'12687510',
                    '12687731',
                    '12687732'}; %Hewitt et al. (2012) 2D
    case 'hewitt_3D'

        slurm_num={'12687744',
                  '12687755',
                  '12687756'}; %Hewitt et al. (2014) 3D
        %'12640113', %%This case is wierd... not right
end

flag.print=0;
flag.visible=0;
flag.video=0;
for slurm_ind=1:length(slurm_num)
    h5_name=['C:\Data\dedalus\dedalus_',...
        slurm_num{slurm_ind},...
            '\analysis\analysis_s1.h5'];
    dedalus_post_my{slurm_ind}=dedalus_post(h5_name,flag);
    dedalus_post_my{slurm_ind}.uvw_hewitt=0;
    dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.dedalus_post_bvp();
    dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.bvp_plot;
    if dedalus_post_my{slurm_ind}.dy_T_mean<0
        background_T=1-dedalus_post_my{slurm_ind}.z_list;
    elseif dedalus_post_my{slurm_ind}.dy_T_mean>0
        background_T=dedalus_post_my{slurm_ind}.z_list;
    end
    data_T{slurm_ind}.y=dedalus_post_my{slurm_ind}.z_list;
    data_T{slurm_ind}.x=dedalus_post_my{slurm_ind}.T_0+background_T;
    data_T{slurm_ind}.x=(data_T{slurm_ind}.x-min(data_T{slurm_ind}.x))/(max(data_T{slurm_ind}.x)-min(data_T{slurm_ind}.x));
    switch group_name
        case {'HB_porous_diffusive_BC','HB_porous_finger_BC'}
          if dedalus_post_my{slurm_ind}.dy_S_mean<0
              background_S=1-dedalus_post_my{slurm_ind}.z_list;
          elseif dedalus_post_my{slurm_ind}.dy_S_mean>0
              background_S=dedalus_post_my{slurm_ind}.z_list;
          end
          data_S{slurm_ind}.y=dedalus_post_my{slurm_ind}.z_list;
          data_S{slurm_ind}.x=dedalus_post_my{slurm_ind}.S_0+background_S;
          data_S{slurm_ind}.x=(data_S{slurm_ind}.x-min(data_S{slurm_ind}.x))/(max(data_S{slurm_ind}.x)-min(data_S{slurm_ind}.x));            
    end
end
switch group_name
    case 'hewitt_2D'
        plot_config.legend_list={1,'$Ra_T=10000$','$Ra_T=20000$','$Ra_T=40000$'};
        plot_config.label_list={1,'$\bar{T}_0+1-z$','$z$'};
        plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'T_0.png'];
        plot_config.fontsize_legend=18;
        plot_config.Markerindex=3;
        plot_config.user_color_style_marker_list={'r-','g--','b-.'};
        plot_config.print_size=[1,600,1000];
        plot_line(data_T,plot_config);
        plot_config.xlim_list=[1,0.45,0.55];
        plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'T_0_local.png'];
        plot_line(data_T,plot_config);
    case 'hewitt_3D'
        plot_config.legend_list={1,'$Ra_T=4000$','$Ra_T=8000$','$Ra_T=16000$'};
        plot_config.label_list={1,'$\bar{T}_0+1-z$','$z$'};
        plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'T_0.png'];
        plot_config.fontsize_legend=18;
        plot_config.Markerindex=3;
        plot_config.user_color_style_marker_list={'r-','g--','b-.'};
        plot_config.print_size=[1,600,1000];
        plot_line(data_T,plot_config);
        plot_config.xlim_list=[1,0.45,0.55];
        plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'T_0_local.png'];
        plot_line(data_T,plot_config);
    case 'HB_porous_thermal_BC'
        plot_config.label_list={1,'$\bar{T}_0+1-z$','$z$'};
        plot_config.Markerindex=3;
        plot_config.print_size=[1,600,1000];
        marker_num=5;
        marker_ind=round(linspace(1,length(data_T{1}.x),marker_num));
        marker_ind=marker_ind(2:end-1);
        for i=2:5
            data_T{i+4}.x=data_T{i}.x(marker_ind);
            data_T{i+4}.y=data_T{i}.y(marker_ind);
        end
        plot_config.user_color_style_marker_list={'k-','b--','b--','r--','r--','bsquare','bo','r^','rv'};
        plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'T_0.png'];
        plot_config.print=1;
        plot_line(data_T,plot_config);
    case {'HB_porous_diffusive_BC','HB_porous_finger_BC'}
        plot_config.Markerindex=3;
        plot_config.print_size=[1,600,1000];
        marker_num=5;
        marker_ind=round(linspace(1,length(data_T{1}.x),marker_num));
        marker_ind=marker_ind(2:end-1);
        for i=2:7
            data_T{i+6}.x=data_T{i}.x(marker_ind);
            data_T{i+6}.y=data_T{i}.y(marker_ind);
            data_S{i+6}.x=data_S{i}.x(marker_ind);
            data_S{i+6}.y=data_S{i}.y(marker_ind);
        end
        plot_config.user_color_style_marker_list={'k-','b--','b--','r--','r--','m--','m--','bsquare','bo','r^','rv','mx','m+'};
        plot_config.print=1;
        plot_config.label_list={1,'$\bar{T}_0+1-z$','$z$'};
        plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'T_0.png'];
        plot_line(data_T,plot_config);
        plot_config.label_list={1,'$\bar{S}_0+1-z$','$z$'};
        plot_config.name=['C:\Figure\DDC_LST\HB_porous_',group_name,'S_0.png'];
        plot_line(data_S,plot_config);
        
        
    case 'HB_porous_finger_BC'
end
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


