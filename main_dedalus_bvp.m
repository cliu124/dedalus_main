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

% group_name='HB_porous_Nu_kx_Ra';
group_name='test';
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
    case 'HB_porous_kx'
        slurm_num={'12761570'}%,...%This is continued from 12639319
                   %'12761797'}; %This is continued from the previous case. 
    case 'HB_porous_Nu_kx_Ra'
        slurm_num={'12784649',...
                    '12784650',...
                    '12784652',...
                    '12784653',...
                    '12784654',...
                    '12784655',...
                    '12784656',...
                    '12784658',...
                    '12784660',...
                    '12784661'};
    case 'hewitt_2D'
        slurm_num={'12687510',
                    '12687731',
                    '12687732'}; %Hewitt et al. (2012) 2D
    case 'wen_chini_2D'
        slurm_num={'12784615'};
        
    case 'hewitt_3D'

        slurm_num={'12687744',
                  '12687755',
                  '12687756'}; %Hewitt et al. (2014) 3D
        %'12640113', %%This case is wierd... not right
    case 'HB_benard_kx'
        slurm_num={'12760884'}; %Ra=10^6, kx=1~20, reproducing figure 7 of Toomre (1977)
    case 'HB_benard_Ra'
        slurm_num={'12761067'}; %kx=1, Ra=10^6~10^10, reproducing figure 4 of Toomre (1977)
    
    case 'herring_1963_free'
        slurm_num={'12761527',...Ra=4*10^3, kx=0.8pi
            '12761528',... Ra=10^4, kx=pi
            '12761529',...Ra=10^5, kx=1.5pi
            '12761536',...Ra=10^6, kx=1.5pi
            '12761538',...Ra=10^6, kx=6pi
            '12761539',...Ra=10^6, kx=9pi 
                    }
    case 'herring_1964_rigid'  
        %slurm_num={'12760848'}; continuation from this case, Ra_T=10^6,
        %kx=1
        slurm_num={'12761520',...Ra=4*10^3, kx=3
                '12761519',...Ra=10^4, kx=3
                '12761518',... Ra=10^5, kx=5
                '12761515'%Ra=10^6, kx=5.25
                }

    case 'test'
        slurm_num={'12804664'};
        %slurm_num={'12760848'}; %Ra=10^6, kx=1
        
                %'12760763',
                %'12760764'};
        %slurm_num={'12760848',%for a range of Ra from 10^4 to 10^6
        %        '12760884'};%for Ra=10^6, but wavenumber from 1 to 20
        
end

flag.print=0;
flag.visible=0;
flag.video=0;
for slurm_ind=1:length(slurm_num)
    content=dir(['C:\Data\dedalus\dedalus_',...
        slurm_num{slurm_ind}]);
    content(1:3)=[];%remove these...
    for content_ind=1:length(content)
        if strcmp(content(content_ind).name(1:8),'analysis')
            h5_name=['C:\Data\dedalus\dedalus_',...
            slurm_num{slurm_ind},...
                '\',content(content_ind).name,'\analysis_s1.h5'];

            dedalus_post_my{slurm_ind,content_ind}=dedalus_post(h5_name,flag);
            dedalus_post_my{slurm_ind,content_ind}.uvw_hewitt=0;
            dedalus_post_my{slurm_ind,content_ind}=dedalus_post_my{slurm_ind,content_ind}.dedalus_post_bvp();
            dedalus_post_my{slurm_ind,content_ind}=dedalus_post_my{slurm_ind,content_ind}.bvp_plot;
            if dedalus_post_my{slurm_ind,content_ind}.dy_T_mean<0
                background_T=1-dedalus_post_my{slurm_ind,content_ind}.z_list;
            elseif dedalus_post_my{slurm_ind,content_ind}.dy_T_mean>0
                background_T=dedalus_post_my{slurm_ind,content_ind}.z_list;
            end
            data_T{slurm_ind,content_ind}.y=dedalus_post_my{slurm_ind,content_ind}.z_list;
            data_T{slurm_ind,content_ind}.x=dedalus_post_my{slurm_ind,content_ind}.T_0+background_T;
            data_T{slurm_ind,content_ind}.x=(data_T{slurm_ind,content_ind}.x-min(data_T{slurm_ind,content_ind}.x))/(max(data_T{slurm_ind,content_ind}.x)-min(data_T{slurm_ind,content_ind}.x));
            switch group_name
                case 'hewitt_2D'
                    data_z{1}.x(slurm_ind)=dedalus_post_my{slurm_ind}.Ra_T;
                    data_z{1}.y(slurm_ind)=dedalus_post_my{slurm_ind}.z_T_BL;
                    data_z{2}.x(slurm_ind)=dedalus_post_my{slurm_ind}.Ra_T;
                    data_z{2}.y(slurm_ind)=dedalus_post_my{slurm_ind}.z_T_rms_max;
                    %data_T_BL{slurm_ind}.x=deda
                    %data_S_BL{slurm_ind}.x=dedalus_post_my{slurm_ind}.Ra_T;
                    %data_S_BL{slurm_ind}.y=dedalus_post_my{slurm_ind}.S_BL;
                    
                case 'hewitt_3D'
                
                case {'HB_porous_diffusive_BC','HB_porous_finger_BC'}
                  if dedalus_post_my{slurm_ind,content_ind}.dy_S_mean<0
                      background_S=1-dedalus_post_my{slurm_ind,content_ind}.z_list;
                  elseif dedalus_post_my{slurm_ind,content_ind}.dy_S_mean>0
                      background_S=dedalus_post_my{slurm_ind,content_ind}.z_list;
                  end
                  data_S{slurm_ind,content_ind}.y=dedalus_post_my{slurm_ind,content_ind}.z_list;
                  data_S{slurm_ind,content_ind}.x=dedalus_post_my{slurm_ind,content_ind}.S_0+background_S;
                  data_S{slurm_ind,content_ind}.x=(data_S{slurm_ind,content_ind}.x-min(data_S{slurm_ind,content_ind}.x))/(max(data_S{slurm_ind,content_ind}.x)-min(data_S{slurm_ind,content_ind}.x));            
                case 'HB_benard_kx'
                    data_Nu{1}.x(content_ind)=dedalus_post_my{slurm_ind,content_ind}.kx;
                    data_Nu{1}.y(content_ind)=dedalus_post_my{slurm_ind,content_ind}.Nu(1);
                case 'HB_benard_Ra'
                    data_Nu{1}.x(content_ind)=dedalus_post_my{slurm_ind,content_ind}.Ra_T;
                    data_Nu{1}.y(content_ind)=dedalus_post_my{slurm_ind,content_ind}.Nu(1);
                case 'HB_porous_kx'
                    data_Nu{1}.x(content_ind)=dedalus_post_my{slurm_ind,content_ind}.kx;
                    data_Nu{1}.y(content_ind)=dedalus_post_my{slurm_ind,content_ind}.Nu(1);
%                 case 'HB_porous_Nu_kx_Ra'
%                     data_Nu{1}.x()=dedalus_post_my{slurm_ind,content_ind}.kx;
%                     data_Nu{1}.y(slurm_ind)=dedalus_post_my{slurm_ind,content_ind}.Ra_T;
%                     data_Nu{1}.z=
            end
        end
    end
end
%%plotting
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
        
        plot_config.loglog=[1,1];
        plot_config.xlim_list=0; 
        plot_config.ylim_list=0
        plot_config.label_list={1,'$Ra_T$','$z_T$'};
        [eta,c0]=scaling(data_z{1}.x,data_z{1}.y);
        plot_line(data_z,plot_config);
    
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
    case 'HB_benard_kx'
        
        Nu_kx=dedalus_post_my{1}.get_Nu_kx_Toomre();
        [data_Nu{1}.x,ind]=sort(data_Nu{1}.x);
        data_Nu{1}.y=data_Nu{1}.y(ind);
        data_Nu{2}.x=Nu_kx(:,1);
        data_Nu{2}.y=Nu_kx(:,2);
        plot_config.Markerindex=3;
        plot_config.user_color_style_marker_list={'b-','k*'};
        plot_config.label_list={1,'$k_x$','Nu'};
        plot_config.print_size=[1,1000,900];
        plot_config.name=['C:\Figure\DDC_LST\',group_name,'Nu_kx_Toomre_1977.png'];
        plot_line(data_Nu,plot_config);
        
        
        clear plot_config;
        T_mean_var=1-dedalus_post_my{1}.z_list;
        T_mean_sign='+1-z';
        plot_config.legend_list{1}=1;
        y_ind=round(linspace(1,length(dedalus_post_my{1}.z_list),5));
        for kx_ind=1:6
        	data_T0{kx_ind}.x=dedalus_post_my{2+5*(kx_ind-1)}.T_0+T_mean_var;
            data_T0{kx_ind}.y=dedalus_post_my{2+5*(kx_ind-1)}.z_list;
            plot_config.legend_list{kx_ind+1}=['$k_x=',num2str(round(dedalus_post_my{2+5*(kx_ind-1)}.kx)),'$'];
            data_T0{kx_ind+6}.x=data_T0{kx_ind}.x(y_ind);
            data_T0{kx_ind+6}.y=data_T0{kx_ind}.y(y_ind);
        end
        plot_config.ylim_list=[1,0,1];
        plot_config.xlim_list=[1,0,1];
        plot_config.print_size=[1,900,900];
        plot_config.loglog=[0,0];
        plot_config.label_list={1,['$\bar{T}_0',T_mean_sign,'$'], '$z$'};
        plot_config.Markerindex=3;
        plot_config.user_color_style_marker_list={'k-','k-','b--','b--','r-.','r-.'...
            ,'ksquare','ko','b^','bv','rx','r+'};
        plot_config.fontsize_legend=20;
        plot_config.name=['C:\Figure\DDC_LST\',group_name,'T_0_kx.png'];
        plot_line(data_T0,plot_config);
        
        
    case 'HB_benard_Ra'
        
        Nu_Ra=dedalus_post_my{1}.get_Nu_Ra_Toomre();
        [data_Nu{1}.x,ind]=sort(data_Nu{1}.x);
        data_Nu{1}.y=data_Nu{1}.y(ind);
        data_Nu{2}.x=Nu_Ra(:,1);
        data_Nu{2}.y=Nu_Ra(:,2)+1;
        plot_config.ylim_list=[1,1,100];
        plot_config.ytick_list=[1,1,10,100];
        plot_config.xtick_list=[1,10^5,10^6,10^7,10^8,10^9,10^10];
        plot_config.Markerindex=3;
        plot_config.loglog=[1,1];
        plot_config.user_color_style_marker_list={'b-','k*'};
        plot_config.label_list={1,'$Ra$','Nu'};
        plot_config.print_size=[1,1000,900];
        plot_config.name=['C:\Figure\DDC_LST\',group_name,'Nu_Ra_Toomre_1977.png'];
        plot_line(data_Nu,plot_config);
        
        clear plot_config;
        T_mean_var=1-dedalus_post_my{1}.z_list;
        T_mean_sign='+1-z';
        plot_config.legend_list{1}=1;
        y_ind=round(linspace(1,length(dedalus_post_my{1}.z_list),5));
        for Ra_ind=1:6
        	data_T0{Ra_ind}.x=dedalus_post_my{1+5*(Ra_ind-1)}.T_0+T_mean_var;
            data_T0{Ra_ind}.y=dedalus_post_my{1+5*(Ra_ind-1)}.z_list;
            plot_config.legend_list{Ra_ind+1}=['$Ra_T=',num2str(round(dedalus_post_my{1+5*(Ra_ind-1)}.Ra_T)),'$'];
            data_T0{Ra_ind+6}.x=data_T0{Ra_ind}.x(y_ind);
            data_T0{Ra_ind+6}.y=data_T0{Ra_ind}.y(y_ind);
        end
        plot_config.ylim_list=[1,0,1];
        plot_config.xlim_list=[1,0,1];
        plot_config.print_size=[1,900,900];
        plot_config.loglog=[0,0];
        plot_config.label_list={1,['$\bar{T}_0',T_mean_sign,'$'], '$z$'};
        plot_config.Markerindex=3;
        plot_config.user_color_style_marker_list={'k-','k-','b--','b--','r-.','r-.'...
            ,'ksquare','ko','b^','bv','rx','r+'};
        plot_config.fontsize_legend=20;
        plot_config.name=['C:\Figure\DDC_LST\',group_name,'T_0_Ra.png'];
        plot_line(data_T0,plot_config);
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


