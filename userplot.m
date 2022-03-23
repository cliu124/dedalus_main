function [p,obj]=userplot(p,wnr) % mod of plotsol for Chebychev-diff setup 
plot_config=p.my.plot_config;
p_old=p;
par=p.u(p.nu+1:end);
obj.dy_T_mean=par(7);
obj.dy_S_mean=par(8);
figure(1)
if p.sw.para>2
    %This is for the Hopf bifurcation
    u=p.hopf.y(1:p.nu,:);
    obj.w_hat=u(1:p.np,:);
    obj.T_hat=u(p.np+1:2*p.np,:);
    obj.S_hat=u(2*p.np+1:3*p.np,:);
    obj.T_0=u(3*p.np+1:4*p.np,:);
    obj.S_0=u(4*p.np+1:5*p.np,:);
    obj.z_list=p.x+p.lx;
    obj.t_list=p.hopf.t;
    fig_num=10;
    figure(fig_num+1)
    mesh(obj.t_list,obj.z_list,obj.w_hat);
    figure(fig_num+2)
    mesh(obj.t_list,obj.z_list,obj.T_hat);
    figure(fig_num+3)
    mesh(obj.t_list,obj.z_list,obj.S_hat);
    figure(fig_num+4)
    mesh(obj.t_list,obj.z_list,obj.T_0+obj.dy_T_mean*obj.z_list);
    figure(fig_num+5)
    mesh(obj.t_list,obj.z_list,obj.S_0+obj.dy_S_mean*obj.z_list);    
else
    %steady bifurcation.
    u=p.u(1:p.nu);
    obj.w_hat=u(1:p.np);
    obj.T_hat=u(p.np+1:2*p.np);
    obj.S_hat=u(2*p.np+1:3*p.np);
    obj.T_0=u(3*p.np+1:4*p.np);
    obj.S_0=u(4*p.np+1:5*p.np);
    obj.z_list=p.x+p.lx;
    fig_num=10;
    figure(fig_num+1)
    plot(obj.w_hat,obj.z_list);
    figure(fig_num+2)
    plot(obj.T_hat,obj.z_list);
    figure(fig_num+3)
    plot(obj.S_hat,obj.z_list);
    figure(fig_num+4)
    plot(obj.T_0+obj.dy_T_mean*obj.z_list,obj.z_list);
    figure(fig_num+5)
    plot(obj.S_0+obj.dy_S_mean*obj.z_list,obj.z_list);
end

if plot_config.visible==1 || plot_config.print==1 || plot_config.post==1
    load([p_old.my.folder_name,'/',plot_config.branch_name,'/',...
        plot_config.point_name,'.mat']);
    
    %load the parameters
    par=p.u(p.nu+1:end); 
%     obj.kx=par(1);
%     obj.ky=par(2);
    obj.Ra_T=par(3);
    obj.Ra_S2T=par(4);
    obj.Le=par(5);
    obj.Pr=par(6);
    obj.dy_T_mean=par(7);
    obj.dy_S_mean=par(8);
    obj.tau=1/obj.Le;
    obj.z_list=p.x+p.lx;

    %read the kx kz computed... 
    [obj.kx,obj.ky]=my_kx(p);
    par(1)=obj.kx;
    par(2)=obj.ky;
    p.u(p.nu+1:end)=par;

    %read chebyshev differentiation matrix. 
    [xt, DM] = chebdif(p.np+2, 2);
    D1=DM(2:end-1,2:end-1,1);
    p.mat.D1=D1/p.lx;
    p.mat.D1_full=DM(:,:,1)/p.lx;
    p.mat.D2=DM(2:end-1,2:end-1,2)/p.lx^2;

    if p.sw.para>2
        %Hopf bifurcation
        u=p.hopf.y(1:p.nu,:);
        obj.w_hat=u(1:p.np,:);
        obj.T_hat=u(p.np+1:2*p.np,:);
        obj.S_hat=u(2*p.np+1:3*p.np,:);
        obj.T_0=u(3*p.np+1:4*p.np,:);
        obj.S_0=u(4*p.np+1:5*p.np,:);
        obj.u_tilde=obj.kx*p.mat.D1*obj.w_hat/(obj.kx^2+obj.ky^2);

    else
        %steady bifurcation
        u=p.u(1:p.nu);  
        obj.w_hat=u(1:p.np);
        obj.T_hat=u(p.np+1:2*p.np);
        obj.S_hat=u(2*p.np+1:3*p.np);
        obj.T_0=u(3*p.np+1:4*p.np);
        obj.S_0=u(4*p.np+1:5*p.np);
        obj.u_tilde=obj.kx*p.mat.D1*obj.w_hat/(obj.kx^2+obj.ky^2);
    end
    
%     R_rho_T2S=obj.Ra_T/obj.Ra_S2T;
%     Ra_S=obj.Ra_S2T/obj.tau;
%     
%     obj.kx=2*pi/(2*14.8211*Ra_S^(-0.2428)/R_rho_T2S^(0.25/2));
%     if obj.ky==0
%         obj.u_tilde=D1*obj.w_hat/obj.kx;
%     else
%        obj.u_tilde=D1*obj.w_hat/2/obj.kx; 
%     end

    obj.no_ylabel=plot_config.no_ylabel;
    
    %make a new directory 
    mkdir([p_old.my.folder_name,'/',plot_config.branch_name,'/graph_',...
        plot_config.point_name]);
    ind_case=strfind(p_old.my.folder_name,'/');
    obj.h5_name=[p_old.my.folder_name,'/',plot_config.branch_name,'/graph_',...
        plot_config.point_name,'/',... folder name up to here
        p_old.my.folder_name(ind_case(end-1)+1:ind_case(end)-1),'_'...%This is the case name, extracted from the folder_name
        strrep(plot_config.branch_name,'/','_'),'_',plot_config.point_name,'_benard___'];

    %Set up the sign for the background temperature and salintiy
    if obj.dy_T_mean==-1
        T_mean_var=1+dy_T_mean*obj.z_list;
        T_mean_sign='+1-z';
        dy_T_mean_sign='-1';
    elseif obj.dy_T_mean==1
        T_mean_var=obj.dy_T_mean*obj.z_list;
        T_mean_sign='+z';
        dy_T_mean_sign='+1';
    end

    if obj.dy_S_mean==-1
        S_mean_var=1+obj.dy_S_mean*obj.z_list;
        S_mean_sign='+1-z';
        dy_S_mean_sign='-1';
    elseif obj.dy_S_mean==1
        S_mean_var=obj.dy_S_mean*obj.z_list;
        S_mean_sign='+z';
        dy_S_mean_sign='+1';
    end

    data{1}.x=obj.T_0(:,1)+T_mean_var;
    data{1}.y=obj.z_list;
%     plot_config.fontsize=20;
    plot_config.label_list={1,['$\bar{T}_0',T_mean_sign,'$'], '$z$'};
    plot_config.print_size=[1,500,900];
    plot_config.name=[obj.h5_name(1:end-3),'_HB_','T_0.png'];
    plot_line(data,plot_config);
    
    if p.sw.para>2
        %plot time average quantity
        data{1}.x=mean(obj.T_0,2)+T_mean_var;
        plot_config.name=[obj.h5_name(1:end-3),'_HB_','T_0_hopf_t_ave.png'];
        plot_line(data,plot_config);
        
        %lift as a contour
        data{1}.z=obj.T_0+T_mean_var*ones(1,p.hopf.tl);
        data{1}.x=obj.t_list;
        plot_config.title_list={1,plot_config.label_list{2}};
        plot_config.label_list{2}='$t/T$';
        plot_config.print_size=[1,1000,900];
        plot_config.colormap='jet';
        plot_config.zlim_list=[1,0,1];
        plot_config.ylim_list=[1,0,1];
        plot_config.name=[obj.h5_name(1:end-3),'_HB_','T_0_hopf.png'];
        plot_contour(data,plot_config);
        plot_config.title_list={0};
    end
    %for computing the Nusselt number, use the full matrix. 
    zero_bc=zeros(1,size(obj.T_0,2));
    p.my.Nu=p.mat.D1_full*[zero_bc;obj.T_0;zero_bc]+obj.dy_T_mean;
    clear data


    data{1}.x=obj.S_0(:,1)+S_mean_var;
    data{1}.y=obj.z_list;
    plot_config.label_list={1,['$\bar{S}_0',S_mean_sign,'$'], '$z$'};
    plot_config.print_size=[1,500,900];
    plot_config.name=[obj.h5_name(1:end-3),'_HB_','S_0.png'];
    plot_line(data,plot_config);
    if p.sw.para>2
        %plot the time average quantity
        data{1}.x=mean(obj.S_0,2)+S_mean_var;
        plot_config.name=[obj.h5_name(1:end-3),'_HB_','S_0_hopf_t_ave.png'];
        plot_line(data,plot_config);
        
        %lift as a contour
        data{1}.z=obj.S_0+S_mean_var*ones(1,p.hopf.tl);
        data{1}.x=obj.t_list;
        plot_config.title_list={1,plot_config.label_list{2}};
        plot_config.label_list{2}='$t/T$';
        plot_config.print_size=[1,1000,900];
        plot_config.colormap='jet';
        plot_config.zlim_list=[1,0,1];
        plot_config.ylim_list=[1,0,1];
        plot_config.name=[obj.h5_name(1:end-3),'_HB_','S_0_hopf.png'];
        plot_contour(data,plot_config);
        plot_config.title_list={0};
    end
    %for computing the Nusselt number, use the full matrix. 
    p.my.Nu_S=p.mat.D1_full*[zero_bc;obj.S_0;zero_bc]+obj.dy_S_mean;

    clear data

    data{1}.x=obj.T_hat(:,1);
    data{1}.y=obj.z_list;
    plot_config.label_list={1,'$\widehat{T}$', '$z$'};
    plot_config.print_size=[1,500,900];
    plot_config.name=[obj.h5_name(1:end-3),'_HB_','T_hat.png'];
    plot_line(data,plot_config);
    if p.sw.para>2
        % plot the time average quantity
        data{1}.x=mean(obj.T_hat,2);
        plot_config.name=[obj.h5_name(1:end-3),'_HB_','T_hat_hopf_t_ave.png'];
        plot_line(data,plot_config);
        
        %lift as the contour
        data{1}.z=obj.T_hat;
        data{1}.x=obj.t_list;
        plot_config.title_list={1,plot_config.label_list{2}};
        plot_config.label_list{2}='$t/T$';
        plot_config.print_size=[1,1000,900];
        plot_config.colormap='bluewhitered';
        plot_config.zlim_list=0;
        plot_config.ylim_list=[1,0,1];
        plot_config.name=[obj.h5_name(1:end-3),'_HB_','T_hat_hopf.png'];
        plot_contour(data,plot_config);
        plot_config.title_list={0};
    end
    clear data


    data{1}.x=obj.S_hat(:,1);
    data{1}.y=obj.z_list;
    plot_config.label_list={1,'$\widehat{S}$', '$z$'};
    plot_config.print_size=[1,500,900];
    plot_config.name=[obj.h5_name(1:end-3),'_HB_','S_hat.png'];
    plot_line(data,plot_config);
    if p.sw.para>2
        % plot the time average quantity
        data{1}.x=mean(obj.S_hat,2);
        plot_config.name=[obj.h5_name(1:end-3),'_HB_','S_hat_hopf_t_ave.png'];
        plot_line(data,plot_config);
        
        %lift as the contour
        data{1}.z=obj.S_hat;
        data{1}.x=obj.t_list;
        plot_config.title_list={1,plot_config.label_list{2}};
        plot_config.label_list{2}='$t/T$';
        plot_config.print_size=[1,1000,900];
        plot_config.colormap='bluewhitered';
        plot_config.zlim_list=0;
        plot_config.ylim_list=[1,0,1];
        plot_config.name=[obj.h5_name(1:end-3),'_HB_','S_hat_hopf.png'];
        plot_contour(data,plot_config);
        plot_config.title_list={0};
    end
    clear data


    data{1}.x=obj.w_hat(:,1);
    data{1}.y=obj.z_list;
    plot_config.label_list={1,'$\widehat{w}$', '$z$'};
    plot_config.print_size=[1,500,900];
    plot_config.name=[obj.h5_name(1:end-3),'_HB_','w_hat.png'];
    plot_line(data,plot_config);
    if p.sw.para>2
        % plot the time average quantity
        data{1}.x=mean(obj.w_hat,2);
        plot_config.name=[obj.h5_name(1:end-3),'_HB_','w_hat_hopf_t_ave.png'];
        plot_line(data,plot_config);
        
        %lift as the contour
        data{1}.z=obj.w_hat;
        data{1}.x=obj.t_list;
        plot_config.title_list={1,plot_config.label_list{2}};
        plot_config.label_list{2}='$t/T$';
        plot_config.print_size=[1,1000,900];
        plot_config.colormap='bluewhitered';
        plot_config.zlim_list=0;
        plot_config.ylim_list=[1,0,1];
        plot_config.name=[obj.h5_name(1:end-3),'_HB_','w_hat_hopf.png'];
        plot_contour(data,plot_config);
        plot_config.title_list={0};
    end
    clear data

    data{1}.x=obj.u_tilde(:,1);
    data{1}.y=obj.z_list;
    plot_config.label_list={1,'$\widetilde{u}$', '$z$'};
    plot_config.print_size=[1,500,900];
    plot_config.name=[obj.h5_name(1:end-3),'_HB_','u_tilde.png'];
    plot_line(data,plot_config);
    if p.sw.para>2
        % plot the time average quantity
        data{1}.x=mean(obj.u_tilde,2);
        plot_config.name=[obj.h5_name(1:end-3),'_HB_','u_tilde_hopf_t_ave.png'];
        plot_line(data,plot_config);
        
        %lift as the contour
        data{1}.z=obj.u_tilde;
        data{1}.x=obj.t_list;
        plot_config.title_list={1,plot_config.label_list{2}};
        plot_config.label_list{2}='$t/T$';
        plot_config.print_size=[1,1000,900];
        plot_config.colormap='bluewhitered';
        plot_config.zlim_list=0;
        plot_config.ylim_list=[1,0,1];
        plot_config.name=[obj.h5_name(1:end-3),'_HB_','u_tilde_hopf.png'];
        plot_contour(data,plot_config);
        plot_config.title_list={0};
    end
    clear data

    %Update 2021/12/06, add streamline...
    
    x=linspace(0,2*pi,1000);
    %z_ind_N=length(obj.z_list)/10;
    %z_ind=round(linspace(1,length(obj.z_list),z_ind_N));
    z_ind=1:length(obj.z_list);
    y=obj.z_list(z_ind);
    [data{2}.x,data{2}.y]=meshgrid(x,y);
    %data{2}.y=obj.z_list;
    data{2}.u=obj.u_tilde(z_ind,1)*real(1i*exp(1i*x));
    data{2}.v=obj.w_hat(z_ind,1)*real(exp(1i*x));
    [data{1}.x,data{1}.y]=meshgrid(x,y);
    data{1}.z=NaN*ones(size(data{1}.x));
%     Du=p.mat.D1_full*[0;obj.u_tilde(z_ind,1);0];
%     data{1}.z=Du(2:end-1)*real(1i*exp(1i*x))...
%         -(obj.kx+obj.ky)*obj.w_hat(z_ind,1)*real(1i*exp(1i*x));
%    
    plot_config.xlim_list=[1,0,2*pi];
    plot_config.xtick_list=[1,0,pi/2,pi,3*pi/2,2*pi];
    plot_config.xticklabels_list={1,'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'};
    plot_config.ylim_list=[1,0,1];
    plot_config.label_list={1,'$x k_x$','$z$'};
    plot_config.streamline=1;
    plot_config.user_color_style_marker_list={'k-','b--'};
    plot_config.panel_num=2;
    plot_config.arrow_ratio=0.8;
    plot_config.linewidth=3;
    plot_config.colorbar=0;
    plot_config.name=[obj.h5_name(1:end-3),'_HB_','streamline.png'];
    plot_contour(data,plot_config);
    if obj.no_ylabel
        plot_config.label_list{3}='';
        plot_config.name=[obj.h5_name(1:end-3),'_HB_','streamline_no_ylabel.png'];
        plot_contour(data,plot_config);
    end
    if p.sw.para>2
        %hopf bifurcation... 
        data{2}.u=mean(obj.u_tilde(z_ind,:),2)*real(1i*exp(1i*x));
        data{2}.v=mean(obj.w_hat(z_ind,:),2)*real(exp(1i*x));
        plot_config.name=[obj.h5_name(1:end-3),'_HB_','streamline_hopf_t_ave.png'];
        plot_contour(data,plot_config);
        plot_config.title_list={0};
    end
    
    clear data;
    plot_config.streamline=0;
    plot_config.panel_num=1;
    %Plot the isocontour of T
    % clear data plot_config
    % plot_config.print=p.my.print;
    % plot_config.visible=p.my.visible;

    x=linspace(0,2*pi,1000);
    z_ind=1:length(obj.z_list);
    y=obj.z_list(z_ind);
    [data{1}.x,data{1}.y]=meshgrid(x,y);
    data{1}.z=T_mean_var+obj.T_0(z_ind,1)+obj.T_hat(z_ind,1)*2*cos(x);
    %data{2}.y=obj.z_list;
    %             data{1}.u=obj.u_tilde(z_ind)*real(1i*exp(1i*x));
    %             data{1}.v=obj.w_hat(z_ind)*real(exp(1i*x));
    plot_config.xlim_list=[1,0,2*pi];
    plot_config.xtick_list=[1,0,pi/2,pi,3*pi/2,2*pi];
    plot_config.xticklabels_list={1,'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'};
    plot_config.ylim_list=[1,0,1];
    plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
    plot_config.label_list={1,'$x k_x$','$z$'};
    plot_config.contour_line=0;
    plot_config.colorbar=1;
    %plot_config.user_color_style_marker_list={'k-','b--'};
    plot_config.name=[obj.h5_name(1:end-3),'_HB_','isocontour_T.png'];
    plot_config.print_size=[1,1000,900];
    %             plot_config.colormap='bluewhitered';
    plot_config.zlim_list=[1,0,1];
    plot_config.ztick_list=[1,0,0.2,0.4,0.6,0.8,1];
    plot_config.colormap='jet';
    plot_contour(data,plot_config);
    if obj.no_ylabel
        plot_config.label_list{3}='';
        plot_config.name=[obj.h5_name(1:end-3),'_HB_','isocontour_T_no_ylabel.png'];
        plot_contour(data,plot_config);
    end
    if p.sw.para>2
        data{1}.z=T_mean_var+mean(obj.T_0(z_ind,:),2)+mean(obj.T_hat(z_ind,:),2)*2*cos(x);
        plot_config.name=[obj.h5_name(1:end-3),'_HB_','isocontour_T_hopf_t_ave.png'];
        plot_contour(data,plot_config);
        plot_config.title_list={0};
    end

    % clear data plot_config
    % plot_config.print=obj.print;
    % plot_config.visible=obj.visible;
    x=linspace(0,2*pi,1000);
    z_ind=1:length(obj.z_list);
    y=obj.z_list(z_ind);
    [data{1}.x,data{1}.y]=meshgrid(x,y);
    data{1}.z=S_mean_var+obj.S_0(z_ind,1)+obj.S_hat(z_ind,1)*2*cos(x);
    plot_config.xlim_list=[1,0,2*pi];
    plot_config.xtick_list=[1,0,pi/2,pi,3*pi/2,2*pi];
    plot_config.xticklabels_list={1,'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'};
    plot_config.ylim_list=[1,0,1];
    plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
    plot_config.label_list={1,'$x k_x$','$z$'};
    plot_config.contour_line=0;
    plot_config.colorbar=1;
    %             plot_config.colormap='bluewhitered';
    plot_config.zlim_list=[1,0,1];
    plot_config.ztick_list=[1,0,0.2,0.4,0.6,0.8,1];
    %plot_config.user_color_style_marker_list={'k-','b--'};
    plot_config.name=[obj.h5_name(1:end-3),'_HB_','isocontour_S.png'];
    plot_config.print_size=[1,1000,900];
    plot_contour(data,plot_config);
    if obj.no_ylabel
        plot_config.label_list{3}='';
        plot_config.name=[obj.h5_name(1:end-3),'_HB_','isocontour_S_no_ylabel.png'];
        plot_contour(data,plot_config);
    end
    if p.sw.para>2
        data{1}.z=S_mean_var+mean(obj.S_0(z_ind,:),2)+mean(obj.S_hat(z_ind,:),2)*2*cos(x);
        plot_config.name=[obj.h5_name(1:end-3),'_HB_','isocontour_S_hopf_t_ave.png'];
        plot_contour(data,plot_config);
        plot_config.title_list={0};
    end
        
    %also output the obj to p.my.....
    p.my.obj=obj;
    
    p_old.my.folder_name=[p_old.my.folder_name,'/',plot_config.branch_name,'/graph_',...
        plot_config.point_name,'/'];
    disp_par(p_old);
    
end





