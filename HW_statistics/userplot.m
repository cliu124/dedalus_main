function [p,obj]=userplot(p,wnr) % mod of plotsol for Chebychev-diff setup 

switch p.my.variable_version
    case {'real_minimal',... only have five  variable, w_hat, T_hat, S_hat, T_0, S_0
            'complex_minimal',... only have 8 variable, add imag of the hat, w_hat, T_hat, S_hat, T_0, S_0, w_hat_imag, T_hat_imag, S_hat_imag
            'complex_zonal',... add one variable for the zonal flow, that should be the 9th variable
            'complex_zonal_with_omega_z'... add two more variable, that should be 11 variable, (most comprehensive one), w_hat, T_hat, S_hat, T_0, S_0, w_hat_imag, T_hat_imag, S_hat_imag, U_0, omega_z_hat, omega_z_hat_imag
            'real_with_omega_z',...A short version to test the hexagon planform, %This should be size variable, w_hat, T_hat, S_hat, T_0, S_0, omega_z_hat
            'complex_with_omega_z'...A short version to test the hexagon planform, but  with the complex one. This should have 10 variable, w_hat, T_hat, S_hat, T_0, S_0, w_hat_imag, T_hat_
            }

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
        obj.U_0=u(8*p.np+1:9*p.np);
        obj.omega_z_hat=u(1+9*p.np:10*p.np);

        obj.z_list=p.x;%+p.lx;
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
        obj.z_list=p.x;%;+p.lx;
        fig_num=10;
        figure(fig_num+1)
        plot(obj.w_hat,obj.z_list);
        figure(fig_num+2)
        plot(obj.T_hat,obj.z_list);
        figure(fig_num+3)
        plot(obj.S_hat,obj.z_list);
        figure(fig_num+4)
        plot(obj.T_0+obj.dy_T_mean*obj.z_list,obj.z_list,'*');
        figure(fig_num+5)
        plot(obj.S_0+obj.dy_S_mean*obj.z_list,obj.z_list,'*');
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
        obj.z_list=p.x;%+p.lx;

        %read the kx kz computed... 
        [obj.kx,obj.ky]=my_kx(p,p.u);
        par(1)=obj.kx;
        par(2)=obj.ky;
        p.u(p.nu+1:end)=par;

        %read chebyshev differentiation matrix. 
        %Update 2022/05/12, modify the differential matrix, keep that for the
        %boundary values.. Need to be
       % error('Update this differential matrix');
        if  strcmp(p.my.z_basis_mode,'Chebyshev')
            if size(p.my.grid,1)==1 
                p.x=[];
                p.mat.D1=[]; p.mat.D2=[]; p.mat.D3=[]; p.mat.D4=[];
                p.mat.Iw=[];
                grid_ind=1;
    %             [xt, DM] = chebdif(p.my.grid(grid_ind,3), 4);%get the chebdif up to fourth order            
                [xt, DM] = chebdif(p.np, 4);%get the chebdif up to fourth order                        
                lx_original=max(xt)-min(xt); %the domain of chebyshev grid, should be 2
                scaling=lx_original/(p.my.grid(grid_ind,2)-p.my.grid(grid_ind,1)); %get scaling factor, 2/(zeta_2-zeta_1)
                for derivative_ind=1:4
                    D_rescale=DM(:,:,derivative_ind)*scaling^derivative_ind;%rescale the differential matrix using the scaling factor
                    %p.mat.(['sub',num2str(grid_ind),'_D',num2str(derivative_ind)])=D_rescale;%set up the differential matrix of each sub-domain, will be used to setup the continuity between each domain
                    p.mat.(['D',num2str(derivative_ind)])=blkdiag(D_rescale,p.mat.(['D',num2str(derivative_ind)]));%set up the whole differential matrix, taking the block diag
                end
                xt_scale=(p.my.grid(grid_ind,2)-p.my.grid(grid_ind,1))*xt/2+(p.my.grid(grid_ind,2)+p.my.grid(grid_ind,1))/2;
                p.x=[xt_scale;p.x];
    %             Iw=[];
                [~,Iw]=clencurt(p.my.grid(grid_ind,3));
                p.mat.Iw=[Iw'/scaling];%setup the integration weight
            else
                error('Not supported!!!');
            end
        elseif strcmp(p.my.z_basis_mode,'Fourier')
             if size(p.my.grid,1)==1
                p.x=[];
                p.mat.D1=[]; p.mat.D2=[]; p.mat.D3=[]; p.mat.D4=[];
                grid_ind=1;
    %             [xt, DM] = chebdif(p.my.grid(grid_ind,3), 4);%get the chebdif up to fourth order            
                [xt, DM(:,:,1)] = fourdif(p.np, 1);%get the chebdif up to fourth order                        
                [~, DM(:,:,2)] = fourdif(p.np, 2);%get the chebdif up to fourth order                        
                [~, DM(:,:,3)] = fourdif(p.np, 3);%get the chebdif up to fourth order                        
                [~, DM(:,:,4)] = fourdif(p.np, 4);%get the chebdif up to fourth order                        

                lx_original=2*pi;%max(xt)-min(xt); %the domain of chebyshev grid, should be 2
                scaling=lx_original/(p.my.grid(grid_ind,2)-p.my.grid(grid_ind,1)); %get scaling factor, 2/(zeta_2-zeta_1)
                for derivative_ind=1:4
                    D_rescale=DM(:,:,derivative_ind)*scaling^derivative_ind;%rescale the differential matrix using the scaling factor
                    %p.mat.(['sub',num2str(grid_ind),'_D',num2str(derivative_ind)])=D_rescale;%set up the differential matrix of each sub-domain, will be used to setup the continuity between each domain
                    p.mat.(['D',num2str(derivative_ind)])=blkdiag(D_rescale,p.mat.(['D',num2str(derivative_ind)]));%set up the whole differential matrix, taking the block diag
                end
                %rescale to thave the same spacing 
                xt_scale=(p.my.grid(grid_ind,2)-p.my.grid(grid_ind,1))*xt/(lx_original);%+(p.my.grid(grid_ind,2)+p.my.grid(grid_ind,1))/2;

                %shift the minimal point in xt_scale to zero and then shift to
                %p.my.grid(grid_ind,1)
                xt_scale=xt_scale-min(xt_scale)+p.my.grid(grid_ind,1);
    %             xt_scale=
                p.x=[xt_scale;p.x];

                dx=diff(p.x); dx=dx(1);
    %             [~,Iw]=clencurt(p.my.grid(grid_ind,3)-1);
                p.mat.Iw=[dx*ones(size(p.x))];%setup the integration weight

                %[~,Iw]=clencurt(p.my.grid(grid_ind,3)-1);
                %p.mat.Iw=[Iw'/scaling;p.mat.Iw];%setup the integration weight
            else
                error('Not supported!!!');

            end

        end


        if p.sw.para>2
            %Hopf bifurcation
            u=p.hopf.y(1:p.nu,:);
            obj.w_hat=u(1:p.np,:);
            obj.T_hat=u(p.np+1:2*p.np,:);
            obj.S_hat=u(2*p.np+1:3*p.np,:);
            obj.T_0=u(3*p.np+1:4*p.np,:);
            obj.S_0=u(4*p.np+1:5*p.np,:);
            obj.U_0=u(8*p.np+1:9*p.np,:);
            obj.omega_z_hat=u(1+9*p.np:10*p.np,:);
            obj.u_tilde=obj.kx*p.mat.D1*obj.w_hat/(obj.kx^2+obj.ky^2);
            plot_complex=0;
        else
            %steady bifurcation
            u=p.u(1:p.nu);  
            obj.w_hat=u(1:p.np);
            obj.T_hat=u(p.np+1:2*p.np);
            obj.S_hat=u(2*p.np+1:3*p.np);
            obj.T_0=u(3*p.np+1:4*p.np);
            obj.S_0=u(4*p.np+1:5*p.np);
            obj.U_0=u(8*p.np+1:9*p.np);
            obj.omega_z_hat=u(1+9*p.np:10*p.np);

            if obj.kx>=0
                obj.w_hat_imag=u(5*p.np+1:6*p.np);
                obj.T_hat_imag=u(6*p.np+1:7*p.np);
                obj.S_hat_imag=u(7*p.np+1:8*p.np);
                obj.omega_z_hat_imag=u(1+10*p.np:11*p.np);
            else
                obj.kx=-obj.kx;
                obj.ky=-obj.ky;
                obj.w_hat_imag=-u(5*p.np+1:6*p.np);
                obj.T_hat_imag=-u(6*p.np+1:7*p.np);
                obj.S_hat_imag=-u(7*p.np+1:8*p.np);
                obj.omega_z_hat_imag=-u(1+10*p.np:11*p.np);
            end

            dw_weight=1;
            obj.u_tilde=dw_weight*obj.kx*p.mat.D1*obj.w_hat/(obj.kx^2+obj.ky^2)...
                    -obj.ky*obj.omega_z_hat/(obj.kx^2+obj.ky^2);
            obj.u_tilde_imag=dw_weight*obj.kx*p.mat.D1*obj.w_hat_imag/(obj.kx^2+obj.ky^2)...
                    -obj.ky*obj.omega_z_hat_imag/(obj.kx^2+obj.ky^2);

            obj.v_tilde=dw_weight*obj.ky*p.mat.D1*obj.w_hat/(obj.kx^2+obj.ky^2)...
                    +obj.kx*obj.omega_z_hat/(obj.kx^2+obj.ky^2);
            obj.v_tilde_imag=dw_weight*obj.ky*p.mat.D1*obj.w_hat_imag/(obj.kx^2+obj.ky^2)...
                    +obj.kx*obj.omega_z_hat_imag/(obj.kx^2+obj.ky^2);
            %phase angle, in 


            %Update 2022/06/01
            %new method to compute phase of harmonics. This will firstly remove
            %any solution that is close to zero... in machine precision as then
            %the phase will be arbitrary and not accurate.
            phase_list={'w_phase','T_phase','S_phase',...
                'u_phase','v_phase'};
            var_list={'w_hat','T_hat','S_hat','u_tilde','v_tilde'}
            for phase_ind=1:length(phase_list)
               ind_non_zero=find(abs(obj.(var_list{phase_ind})+1i*obj.([var_list{phase_ind},'_imag']))>1e-8);
               obj.(phase_list{phase_ind})=zeros(length(obj.z_list),1);
               obj.(phase_list{phase_ind})(ind_non_zero)=...
                   atan(obj.([var_list{phase_ind},'_imag'])(ind_non_zero)./...
                   obj.(var_list{phase_ind})(ind_non_zero));
            end
    %         obj.w_phase=atan(obj.w_hat_imag./obj.w_hat);
    %         obj.T_phase=atan(obj.T_hat_imag./obj.T_hat);
    %         obj.S_phase=atan(obj.S_hat_imag./obj.S_hat);
    %         obj.u_phase=atan(obj.u_tilde_imag./obj.u_tilde);
    %         obj.v_phase=atan(obj.v_tilde_imag./obj.v_tilde);

            %Update 2022/05/12, only concern about the phase in between
            if sum(abs(obj.w_phase(2:end-1)-obj.T_phase(2:end-1)))<0.001 & sum(abs(obj.w_phase(2:end-1)-obj.S_phase(2:end-1)))<0.001
                %all of these phase are the same... just translate to the zero
                %phase, and make w_hat, T_hat, and S_hat as the amplitude of
                %these variable, but also keep the sign information.
                plot_complex=0;
                obj.w_hat=obj.w_hat./cos(obj.w_phase);
                obj.T_hat=obj.T_hat./cos(obj.T_phase);
                obj.S_hat=obj.S_hat./cos(obj.S_phase);
                obj.u_tilde=obj.u_tilde./cos(obj.u_phase);
                obj.v_tilde=obj.v_tilde./cos(obj.v_phase);

                %replace any NaN as zero. That is due to the division by zero
                %error leading to NaN
                %Update 2022/05/12
                obj.w_hat(isnan(obj.w_hat))=0;
                obj.T_hat(isnan(obj.T_hat))=0;
                obj.S_hat(isnan(obj.S_hat))=0;
                obj.u_tilde(isnan(obj.u_tilde))=0;
                obj.v_tilde(isnan(obj.v_tilde))=0;
            else
                %The phase at each y location are not the same.. just add them
                %together.
                plot_complex=1;
                obj.w_hat=obj.w_hat+1i*obj.w_hat_imag;
                obj.T_hat=obj.T_hat+1i*obj.T_hat_imag;
                obj.S_hat=obj.S_hat+1i*obj.S_hat_imag;
                obj.u_tilde=obj.u_tilde+1i*obj.u_tilde_imag;
                obj.v_tilde=obj.v_tilde+1i*obj.v_tilde_imag;
                obj.omega_z_hat=obj.omega_z_hat+1i*obj.omega_z_hat_imag;

            end

        end

    %     if strcmp(p.my.z_basis_mode,'Fourier')
    %        %also modify 
    %        obj.z_list=[obj.z_list;obj.z_list(end)+obj.z_list(2)];
    %        var_name_list={'w_hat','T_hat','S_hat','T_0','S_0','U_0','omega_z_hat','u_tilde','v_tilde'};
    %        for var_name_ind=1:length(var_name_list)
    %            var_name=var_name_list{var_name_ind};
    %            obj.(var_name)=[obj.(var_name);obj.(var_name)(1)];
    %        end
    %     end


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
            T_mean_var=1+obj.dy_T_mean*obj.z_list;
            T_mean_sign='1-z+';
            dy_T_mean_sign='-1';
        elseif obj.dy_T_mean==1
            T_mean_var=obj.dy_T_mean*obj.z_list;
            T_mean_sign='z+';
            dy_T_mean_sign='+1';
        elseif obj.dy_T_mean==0
            T_mean_var=0;
            T_mean_sign='';
            dy_T_mean_sign='';
        end

        if p.my.flux_T
            dy_T_mean_q=p.u(p.nu);
            obj.dy_T_mean=dy_T_mean_q;
            T_mean_var=1+dy_T_mean_q*obj.z_list;
            T_mean_sign='1+\bar{\mathcal{T}}_{z,q}z+';
        end
        if obj.dy_S_mean==-1
            S_mean_var=1+obj.dy_S_mean*obj.z_list;
            S_mean_sign='1-z+';
            dy_S_mean_sign='-1';
        elseif obj.dy_S_mean==1
            S_mean_var=obj.dy_S_mean*obj.z_list;
            S_mean_sign='z+';
            dy_S_mean_sign='+1';
        elseif obj.dy_S_mean==0
            S_mean_var=0;
            S_mean_sign='';
            dy_S_mean_sign='';
        end

        data{1}.x=obj.T_0(:,1)+T_mean_var;
        data{1}.y=obj.z_list;
    %     plot_config.fontsize=20;
        plot_config.linewidth=3;
        plot_config.label_list={1,['$',T_mean_sign,'\bar{T}_0','$'], '$z$'};
        plot_config.print_size=[1,500,900];
        plot_config.name=[obj.h5_name(1:end-3),'_HB_','T_0.png'];
        if (max(obj.z_list)-min(obj.z_list))>1
            plot_config.ytick_list=0;
        else
            plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
        end
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
            plot_config.xlim_list=[1,0,1];
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','T_0_hopf.png'];
            plot_contour(data,plot_config);
            plot_config.title_list={0};
        end
        %for computing the Nusselt number, use the full matrix. 
        zero_bc=zeros(1,size(obj.T_0,2));
        p.my.Nu=p.mat.D1*[obj.T_0]+obj.dy_T_mean;
        clear data


        data{1}.x=obj.S_0(:,1)+S_mean_var;
        data{1}.y=obj.z_list;
        plot_config.label_list={1,['$',S_mean_sign,'\bar{S}_0','$'], '$z$'};
        plot_config.print_size=[1,500,900];
        plot_config.name=[obj.h5_name(1:end-3),'_HB_','S_0.png'];
        if (max(obj.z_list)-min(obj.z_list))>1
            plot_config.ytick_list=0;
        else
            plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
        end
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
            plot_config.xlim_list=[1,0,1];
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','S_0_hopf.png'];
            plot_contour(data,plot_config);
            plot_config.title_list={0};
        end
        %for computing the Nusselt number, use the full matrix. 
        p.my.Nu_S=p.mat.D1*[obj.S_0]+obj.dy_S_mean;

        %%Update 2022/04/04, add the plotting of the large scale shear (zonal flow)
        data{1}.x=obj.U_0(:,1);
        data{1}.y=obj.z_list;
    %     plot_config.fontsize=20;
        plot_config.label_list={1,['$\bar{U}_0$'], '$z$'};
        plot_config.print_size=[1,500,900];
        plot_config.name=[obj.h5_name(1:end-3),'_HB_','U_0.png'];
        if (max(obj.z_list)-min(obj.z_list))>1
            plot_config.ytick_list=0;
        else
            plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
        end
       plot_line(data,plot_config);
        if p.sw.para>2

            %lift as the contour
            data{1}.z=obj.U_0;
            data{1}.x=obj.t_list;
            plot_config.title_list={1,plot_config.label_list{2}};
            plot_config.label_list{2}='$t/T$';
            plot_config.print_size=[1,1000,900];
            plot_config.colormap='bluewhitered';
            plot_config.zlim_list=0;
            plot_config.ylim_list=[1,0,1];
            plot_config.xlim_list=[1,0,1];
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','U_0_hopf.png'];
            plot_contour(data,plot_config);
            plot_config.title_list={0};
        end

        clear data

        data{1}.x=real(obj.T_hat(:,1));
        data{1}.y=obj.z_list;
        plot_config.label_list={1,'$\widehat{T}$', '$z$'};
        plot_config.print_size=[1,500,900];
        plot_config.name=[obj.h5_name(1:end-3),'_HB_','T_hat.png'];
        if plot_complex
            data{2}.x=imag(obj.T_hat(:,1));
            data{2}.y=obj.z_list;
            plot_config.legend_list={1,'$\mathcal{R}e[\cdot]$','$\mathcal{I}m[\cdot]$'};
            plot_config.fontsize_legend=18;
        end
        if (max(obj.z_list)-min(obj.z_list))>1
            plot_config.ytick_list=0;
        else
            plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
        end
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
            plot_config.xlim_list=[1,0,1];
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','T_hat_hopf.png'];
            plot_contour(data,plot_config);
            plot_config.title_list={0};
        end
        clear data


        data{1}.x=real(obj.S_hat(:,1));
        data{1}.y=obj.z_list;
        plot_config.label_list={1,'$\widehat{S}$', '$z$'};
        plot_config.print_size=[1,500,900];
        plot_config.name=[obj.h5_name(1:end-3),'_HB_','S_hat.png'];
        if plot_complex
            data{2}.x=imag(obj.S_hat(:,1));
            data{2}.y=obj.z_list;
            plot_config.legend_list={1,'$\mathcal{R}e[\cdot]$','$\mathcal{I}m[\cdot]$'};
            plot_config.fontsize_legend=18;
        end
        x_max=max([abs(real(obj.S_hat(:,1)));abs(imag(obj.S_hat(:,1)))]);
        plot_config.xlim_list=[1,-1.1*x_max,1.1*x_max];
        if (max(obj.z_list)-min(obj.z_list))>1
            plot_config.ytick_list=0;
        else
            plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
        end
       plot_line(data,plot_config);
        if obj.no_ylabel
            plot_config.label_list{3}='';
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','S_hat_no_ylabel.png'];
            plot_line(data,plot_config);
        end


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
            plot_config.xlim_list=[1,0,1];
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','S_hat_hopf.png'];
            plot_contour(data,plot_config);
            plot_config.title_list={0};
        end
        clear data


        data{1}.x=real(obj.w_hat(:,1));
        data{1}.y=obj.z_list;
        plot_config.label_list={1,'$\widehat{w}$', '$z$'};
        plot_config.print_size=[1,500,900];
        plot_config.name=[obj.h5_name(1:end-3),'_HB_','w_hat.png'];
        if plot_complex
            data{2}.x=imag(obj.w_hat(:,1));
            data{2}.y=obj.z_list;
            plot_config.legend_list={1,'$\mathcal{R}e[\cdot]$','$\mathcal{I}m[\cdot]$'};
            plot_config.fontsize_legend=18;
        end
        x_max=max([abs(real(obj.w_hat(:,1)));abs(imag(obj.w_hat(:,1)))]);
        plot_config.xlim_list=[1,-1.1*x_max,1.1*x_max];
        %plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
        if (max(obj.z_list)-min(obj.z_list))>1
            plot_config.ytick_list=0;
        else
            plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
        end

        plot_line(data,plot_config);
        if obj.no_ylabel %plot the version without y label for paper writing
            plot_config.label_list{3}='';
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','w_hat_no_ylabel.png'];
            plot_line(data,plot_config);
        end

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
            plot_config.xlim_list=[1,0,1];
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','w_hat_hopf.png'];
            plot_contour(data,plot_config);
            plot_config.title_list={0};
        end
        clear data

        data{1}.x=real(obj.u_tilde(:,1));
        data{1}.y=obj.z_list;
        plot_config.label_list={1,'$\widetilde{u}$', '$z$'};
        plot_config.print_size=[1,500,900];
        plot_config.name=[obj.h5_name(1:end-3),'_HB_','u_tilde.png'];
        if plot_complex
            data{2}.x=imag(obj.u_tilde(:,1));
            data{2}.y=obj.z_list;
            plot_config.legend_list={1,'$\mathcal{R}e[\cdot]$','$\mathcal{I}m[\cdot]$'};
            plot_config.fontsize_legend=18;
        end
        x_max=max([abs(real(obj.u_tilde(:,1)));abs(imag(obj.u_tilde(:,1)))]);
        plot_config.xlim_list=[1,-1.1*x_max,1.1*x_max];
        %plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
        if (max(obj.z_list)-min(obj.z_list))>1
            plot_config.ytick_list=0;
        else
            plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
        end


        plot_line(data,plot_config);
        if obj.no_ylabel %plot the version without y label for paper writing
            plot_config.label_list{3}='';
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','u_tilde_no_ylabel.png'];
            plot_line(data,plot_config);
        end

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

        %Update 2022/04/11, add the plotting of vertical vorticity
        data{1}.x=obj.omega_z_hat(:,1);
        data{1}.y=obj.z_list;
        plot_config.label_list={1,'$\widehat{\zeta}$', '$z$'};
        plot_config.print_size=[1,500,900];
        plot_config.name=[obj.h5_name(1:end-3),'_HB_','omega_z_hat.png'];
        if plot_complex
            data{2}.x=imag(obj.omega_z_hat(:,1));
            data{2}.y=obj.z_list;
            plot_config.legend_list={1,'$\mathcal{R}e[\cdot]$','$\mathcal{I}m[\cdot]$'};
            plot_config.fontsize_legend=18;
        end
        %x_max=max([real(obj.omega_z_hat(:,1));imag(obj.omega_z_hat(:,1))]);
        plot_config.xlim_list=[1,-1.1*x_max,1.1*x_max];
        plot_config.xlim_list=0; plot_config.xtick_list=0;
        %plot_config.xtick_list=[1,-0.02,0,0.02];
        plot_config.xtick_list=[1,-0.1,0,0.1];
%plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
        if (max(obj.z_list)-min(obj.z_list))>1
            plot_config.ytick_list=0;
        else
            plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
        end

        plot_line(data,plot_config);



        %Update 2021/12/06, add streamline...

        x=linspace(0,2*pi,1000);
        %z_ind_N=length(obj.z_list)/10;
        %z_ind=round(linspace(1,length(obj.z_list),z_ind_N));
        z_ind=1:length(obj.z_list);
        y=obj.z_list(z_ind);


        psi_0=zeros(size(obj.U_0));
        for psi_ind=1:length(psi_0)
           psi_0(psi_ind,1)=sum(obj.U_0(psi_ind:end).*p.mat.Iw(psi_ind:end)); 
        end
        %[data{2}.x,data{2}.y]=meshgrid(x,y);
        data{2}.x=x; data{2}.y=y;
        %data{2}.y=obj.z_list;

        %Update 2022/06/02, add the streamfunction corresponding to large-scale
        %shear
        data{2}.z=-psi_0*ones(1,length(x))+2*real(obj.w_hat(z_ind,1)/(1i*obj.kx)*exp(1i*x));

        %old version
    %     data{2}.z=2*real(obj.w_hat(z_ind,1)/(1i*obj.kx)*exp(1i*x));

        %data{2}.u=2*real(obj.u_tilde(z_ind,1)*1i*exp(1i*x));%obj.U_0+
        %data{2}.v=2*real(obj.w_hat(z_ind,1)*exp(1i*x));
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
        plot_config.streamline=2;
        plot_config.user_color_style_marker_list={'k-','r--'};
        plot_config.panel_num=2;
        plot_config.arrow_ratio=0.8;
        plot_config.linewidth=3;
        plot_config.colorbar=0;
%         plot_config.visible=1;
        plot_config.fontsize=28;
        plot_config.print_size=[1,500,900];
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
        plot_config.fontsize=40;

        clear data;
        if obj.kx~=0 && obj.ky~=0
            try obj.z_slices_list=plot_config.z_slices_list;
            catch
                obj.z_slices_list=[0.25,0.5,0.75];
            end
            %plot the velocity vector in the horizontal direction 
            %with vertical vorticity as  the contour... This will be done in a
            %selected wall-normal location.
            for z_slices_ind=1:length(obj.z_slices_list)
                z_slices=obj.z_slices_list(z_slices_ind);

                x=linspace(0,2*pi,20);
                [~,z_ind]=min(abs(z_slices-obj.z_list));
                y=linspace(0,2*pi,18);
                [data{2}.x,data{2}.y]=meshgrid(x,y);
                %get the superposition of the cos(kx x+ky y)+cos(kx x-ky y).
                %This will give a squared planform instead of the inclined roll.
                %note that for v and omega_z, we also need to change the
                %sign...
                data{2}.u=obj.U_0(z_ind,1)...
                    +real(obj.u_tilde(z_ind,1)*1i*exp(1i*y')*exp(1i*x))...
                    +real(obj.u_tilde(z_ind,1)*1i*exp(-1i*y')*exp(1i*x));
                data{2}.v=real(obj.v_tilde(z_ind,1)*1i*exp(1i*y')*exp(1i*x))...
                    -real(obj.v_tilde(z_ind,1)*1i*exp(-1i*y')*exp(1i*x));
                [data{1}.x,data{1}.y]=meshgrid(x,y);
                if max(abs(obj.omega_z_hat(z_ind,1)))<0.01
                    data{1}.z=0*exp(1i*y')*exp(1i*x);
                    plot_config.zlim_list=[1,-0.01,0.01];
                else
                    data{1}.z=real(obj.omega_z_hat(z_ind,1)*exp(1i*y')*exp(1i*x))...
                        -real(obj.omega_z_hat(z_ind,1)*exp(-1i*y')*exp(1i*x));
                    plot_config.zlim_list=0;
                end
                %obj.d_w_hat=p.mat.D1*obj.w_hat;
                %data{1}.z=real(obj.d_w_hat(z_ind,1)*exp(1i*y')*exp(1i*x))...
                %    +real(obj.d_w_hat(z_ind,1)*exp(-1i*y')*exp(1i*x));

                plot_config.xlim_list=[1,0,2*pi];
                plot_config.xtick_list=[1,0,pi/2,pi,3*pi/2,2*pi];
                plot_config.xticklabels_list={1,'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'};
                plot_config.ytick_list=[1,0,pi/2,pi,3*pi/2,2*pi];
                plot_config.yticklabels_list={1,'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'};
                plot_config.ylim_list=[1,0,2*pi];
                plot_config.label_list={1,'$x k_x$','$y k_y$'};
                plot_config.streamline=0;
                plot_config.user_color_style_marker_list={'k-','r--'};
                plot_config.panel_num=2;
                plot_config.arrow_ratio=0.8;
                plot_config.linewidth=3;
                plot_config.colorbar=1;
                plot_config.colormap='bluewhitered';
                plot_config.print_size=[1,1100,900];
                plot_config.name=[obj.h5_name(1:end-3),'_HB_','z_slices=',num2str(z_slices),'.png'];
                plot_contour(data,plot_config);
                if obj.no_ylabel
                    plot_config.label_list{3}='';
                    plot_config.name=[obj.h5_name(1:end-3),'_HB_','z_slices',num2str(z_slices),'_no_ylabel.png'];
                    plot_contour(data,plot_config);
                end
                if p.sw.para>2
                    error('z slices for hopf bifurcation needs to be finished');
        %             %hopf bifurcation... 
        %             data{2}.u=mean(obj.u_tilde(z_ind,:),2)*real(1i*exp(1i*x));
        %             data{2}.v=mean(obj.w_hat(z_ind,:),2)*real(exp(1i*x));
        %             plot_config.name=[obj.h5_name(1:end-3),'_HB_','streamline_hopf_t_ave.png'];
        %             plot_contour(data,plot_config);
        %             plot_config.title_list={0};
                end
            end
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
        data{1}.z=T_mean_var+obj.T_0(z_ind,1)+2*real(obj.T_hat(z_ind,1)*exp(1i*x));
        %data{2}.y=obj.z_list;
        %             data{1}.u=obj.u_tilde(z_ind)*real(1i*exp(1i*x));
        %             data{1}.v=obj.w_hat(z_ind)*real(exp(1i*x));
        plot_config.xlim_list=[1,0,2*pi];
        plot_config.xtick_list=[1,0,pi/2,pi,3*pi/2,2*pi];
        plot_config.xticklabels_list={1,'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'};
        plot_config.ylim_list=[1,0,1];
        %plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
        if (max(obj.z_list)-min(obj.z_list))>1
            plot_config.ytick_list=0;
        else
            plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
        end

        plot_config.yticklabels_list={0};
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
        data{1}.z=S_mean_var+obj.S_0(z_ind,1)+2*real(obj.S_hat(z_ind,1)*exp(1i*x));
        plot_config.xlim_list=[1,0,2*pi];
        plot_config.xtick_list=[1,0,pi/2,pi,3*pi/2,2*pi];
        plot_config.xticklabels_list={1,'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'};
        plot_config.ylim_list=[1,0,1];
        %plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
        if (max(obj.z_list)-min(obj.z_list))>1
            plot_config.ytick_list=0;
        else
            plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
        end

        plot_config.yticklabels_list={0};
        plot_config.label_list={1,'$x k_x$','$z$'};
        plot_config.contour_line=0;
        plot_config.colorbar=1;
        %             plot_config.colormap='bluewhitered';
        if strcmp(p.my.z_basis_mode,'Fourier')
            plot_config.zlim_list=0;
            plot_config.ztick_list=0;
        else
            plot_config.zlim_list=[1,0,1];
            plot_config.ztick_list=[1,0,0.2,0.4,0.6,0.8,1];
        end
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
        p.my.folder_name=p_old.my.folder_name;
        disp_par(p);
    end

    %update the 2D configuration here
    case  {'2D_psi_T_S','2D_psi_T_S_decomposed'}
        plot_config=p.my.plot_config;
        p_old=p;
%         
%         figure(101);
%         
%         u_vel_2D=reshape(u_vel,p.nz,p.nx);
%         mesh(u_vel_2D);
%         figure(102);
%         plot(w_vel);
%         figure(103);
%         plot(mean(u_vel_2D,2));
        
        if plot_config.visible==1 || plot_config.print==1 || plot_config.post==1
            load([p_old.my.folder_name,'/',plot_config.branch_name,'/',...
                plot_config.point_name,'.mat']);

            
            p.z=[];
            p.mat.D1=[]; p.mat.D2=[]; p.mat.D3=[]; p.mat.D4=[];
            grid_ind=1;
%             [xt, DM] = chebdif(p.my.grid(grid_ind,3), 4);%get the chebdif up to fourth order            
            [xt, DM(:,:,1)] = fourdif(p.nz, 1);%get the chebdif up to fourth order                        
            %[~, DM(:,:,2)] = fourdif(p.nz, 2);%get the chebdif up to fourth order                        
            %[~, DM(:,:,3)] = fourdif(p.nz, 3);%get the chebdif up to fourth order                        
            %[~, DM(:,:,4)] = fourdif(p.nz, 4);%get the chebdif up to fourth order                        

            lx_original=2*pi;%max(xt)-min(xt); %the domain of chebyshev grid, should be 2
            scaling=lx_original/(p.my.grid(grid_ind,2)-p.my.grid(grid_ind,1)); %get scaling factor, 2/(zeta_2-zeta_1)
            for derivative_ind=1
                D_rescale=DM(:,:,derivative_ind)*scaling^derivative_ind;%rescale the differential matrix using the scaling factor
                %p.mat.(['sub',num2str(grid_ind),'_D',num2str(derivative_ind)])=D_rescale;%set up the differential matrix of each sub-domain, will be used to setup the continuity between each domain
                p.mat.(['D',num2str(derivative_ind)])=blkdiag(D_rescale,p.mat.(['D',num2str(derivative_ind)]));%set up the whole differential matrix, taking the block diag
            end
            %rescale to thave the same spacing 
            xt_scale=(p.my.grid(grid_ind,2)-p.my.grid(grid_ind,1))*xt/(lx_original);%+(p.my.grid(grid_ind,2)+p.my.grid(grid_ind,1))/2;

            %shift the minimal point in xt_scale to zero and then shift to
            %p.my.grid(grid_ind,1)
            xt_scale=xt_scale-min(xt_scale)+p.my.grid(grid_ind,1);
%             xt_scale=
            p.z=[xt_scale;p.z];

%             %add the integration weighting matrix. just the uniform scaled
%             %by dx for fourier basis
%             dz=diff(p.z);
%             dz=dz(1);
%             p.mat.Iw_z=dz*eye(p.nz,p.nz);

            %for 2D, product with kroneck of identify
            for derivative_ind=1
                p.mat.(['D',num2str(derivative_ind),'_z_2D'])=...
                  kron(eye(p.nx,p.nx),p.mat.(['D',num2str(derivative_ind)]));
            end

            %add the horizontal differential matrix, but in [0,2pi]
            %domain
            [p.x, p.mat.D1_2pi] = fourdif(p.nx, 1);%get the chebdif up to fourth order                        
%             [~, p.mat.D2_2pi] = fourdif(p.nx, 2);%get the chebdif up to fourth order                        
%             [~, p.mat.D3_2pi] = fourdif(p.nx, 3);%get the chebdif up to fourth order                        
%             [~, p.mat.D4_2pi] = fourdif(p.nx, 4);%get the chebdif up to fourth order                        
%             dx=diff(p.x);
%             dx=dx(1);
%             p.mat.Iw_x=dx*eye(p.nx,p.nx);
              for derivative_ind=1
                  p.mat.(['D',num2str(derivative_ind),'_2pi_2D'])=...
                    kron(p.mat.(['D',num2str(derivative_ind),'_2pi']),eye(p.nz,p.nz));
              end
            
              %make the matrix that construct the horizontal mean or
              %project the horizontal mean into 2D results
              
            
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
           obj.no_ylabel=plot_config.no_ylabel;

            %obj.z=p.z;
             %%
    %          1
            obj.psi=p.u(1:p.np);%reshape(p.u(1:p.np),p.nz,p.nx);
            obj.T=p.u(1+p.np:2*p.np);%reshape(p.u(1+p.np:2*p.np),p.nz,p.nx);
            obj.S=p.u(1+2*p.np:3*p.np);%reshape(p.u(1+2*p.np:3*p.np),p.nz,p.nx);
            D1z=p.mat.D1_z_2D;
            par=p.u(1+p.nu:end);
            
            obj.kx=par(1);
            D1x=p.mat.D1_2pi_2D*abs(obj.kx);
            obj.z_list=p.z;%+p.lx;
            obj.x_list=p.x/abs(obj.kx);
            obj.w_vel=D1x*obj.psi;
            obj.u_vel=-D1z*obj.psi;

            obj.u_2D=reshape(obj.u_vel,p.nz,p.nx);
            obj.w_2D=reshape(obj.w_vel,p.nz,p.nx);

            obj.T_2D=reshape(obj.T,p.nz,p.nx);
            obj.S_2D=reshape(obj.S,p.nz,p.nx);

            obj.psi_2D=reshape(obj.psi,p.nz,p.nx);

            obj.U_0=mean(obj.u_2D,2);
            obj.T_0=mean(obj.T_2D,2);
            obj.S_0=mean(obj.S_2D,2);
            
            if  strcmp(p.my.variable_version,'2D_psi_T_S_decomposed')
               
                Iz=eye(p.nz,p.nz);
                p.mat.h_mean=1/p.nx*kron(ones(1,p.nx),Iz);
                p.mat.h_mean_2_xz=kron(ones(p.nx,1),Iz);

                obj.psi_0=p.u(1+3*p.np:3*p.np+p.nz);
                obj.U_0=-p.mat.D1*obj.psi_0;
                obj.T_0=p.u(1+3*p.np+p.nz:3*p.np+2*p.nz);
                obj.S_0=p.u(1+3*p.np+2*p.nz:3*p.np+3*p.nz);
                
                obj.psi_p=p.u(1:p.np);
                obj.T_p=p.u(1+p.np:2*p.np);
                obj.S_p=p.u(1+2*p.np:3*p.np);
                
                obj.w_vel=D1x*obj.psi_p;
                obj.u_vel=-D1z*obj.psi_p;
                
                obj.u_2D=reshape(obj.u_vel+p.mat.h_mean_2_xz*obj.U_0,p.nz,p.nx);%add the horizontal mean back
                obj.w_2D=reshape(obj.w_vel,p.nz,p.nx);
                obj.T_2D=reshape(obj.T_p+p.mat.h_mean_2_xz*obj.T_0,p.nz,p.nx);%add horizontal mean back
                obj.S_2D=reshape(obj.S_p+p.mat.h_mean_2_xz*obj.S_0,p.nz,p.nx);%add horizontal mean back
                obj.psi_2D=reshape(obj.psi_p+p.mat.h_mean_2_xz*obj.psi_0,p.nz,p.nx);%add the horizontal mean back
  
            end
            
            
            
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
                T_mean_sign='1-z+';
                dy_T_mean_sign='-1';
            elseif obj.dy_T_mean==1
                T_mean_var=obj.dy_T_mean*obj.z_list;
                T_mean_sign='z+';
                dy_T_mean_sign='+1';
            end

            if obj.dy_S_mean==-1
                S_mean_var=1+obj.dy_S_mean*obj.z_list;
                S_mean_sign='1-z+';
                dy_S_mean_sign='-1';
            elseif obj.dy_S_mean==1
                S_mean_var=obj.dy_S_mean*obj.z_list;
                S_mean_sign='z+';
                dy_S_mean_sign='+1';
            end

            data{1}.x=obj.T_0(:,1)+T_mean_var;
            data{1}.y=obj.z_list;
        %     plot_config.fontsize=20;
            plot_config.label_list={1,['$',T_mean_sign,'\langle T\rangle_h','$'], '$z$'};
            plot_config.print_size=[1,500,900];
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','T_0.png'];
            if (max(obj.z_list)-min(obj.z_list))>1
                plot_config.ytick_list=0;
            else
                plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
            end
            plot_line(data,plot_config);



            data{1}.x=obj.S_0(:,1)+S_mean_var;
            data{1}.y=obj.z_list;
            plot_config.label_list={1,['$',S_mean_sign,'\langle S \rangle_h','$'], '$z$'};
            plot_config.print_size=[1,500,900];
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','S_0.png'];
            plot_config.linewidth=3;
            if (max(obj.z_list)-min(obj.z_list))>1
                plot_config.ytick_list=0;
            else
                plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
            end
            plot_line(data,plot_config);
        
        
            data{1}.x=obj.U_0(:,1);
            data{1}.y=obj.z_list;
        %     plot_config.fontsize=20;
            plot_config.label_list={1,['$\langle u\rangle_h$'], '$z$'};
            plot_config.print_size=[1,500,900];
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','U_0.png'];
            if (max(obj.z_list)-min(obj.z_list))>1
                plot_config.ytick_list=0;
            else
                plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
            end
           plot_line(data,plot_config);

           
           
            clear data;
            plot_config.streamline=0;
            plot_config.panel_num=1;
            %Plot the isocontour of T
            % clear data plot_config
            % plot_config.print=p.my.print;
            % plot_config.visible=p.my.visible;

            %x=linspace(0,2*pi,1000);
            x=obj.x_list;
            z_ind=1:length(obj.z_list);
            y=obj.z_list(z_ind);
            [data{1}.x,data{1}.y]=meshgrid(x,y);
            data{1}.z=T_mean_var*ones(1,p.nx)+obj.T_2D;
            %data{2}.y=obj.z_list;
            %             data{1}.u=obj.u_tilde(z_ind)*real(1i*exp(1i*x));
            %             data{1}.v=obj.w_hat(z_ind)*real(exp(1i*x));
            %plot_config.xlim_list=[1,0,2*pi];
            plot_config.xlim_list=[1,0,max(obj.x_list)];
            %plot_config.xtick_list=[1,0,pi/2,pi,3*pi/2,2*pi];
            %plot_config.xticklabels_list={1,'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'};
            plot_config.ylim_list=[1,0,1];
            %plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
            if (max(obj.z_list)-min(obj.z_list))>1
                plot_config.ytick_list=0;
            else
                plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
            end

            plot_config.yticklabels_list={0};
            plot_config.label_list={1,'$x$','$z$'};
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
            

            % clear data plot_config
            % plot_config.print=obj.print;
            % plot_config.visible=obj.visible;
%             x=linspace(0,2*pi,1000);
            x=obj.x_list
            z_ind=1:length(obj.z_list);
            y=obj.z_list(z_ind);
            [data{1}.x,data{1}.y]=meshgrid(x,y);
            data{1}.z=S_mean_var*ones(1,p.nx)+obj.S_2D;%+obj.S_0(z_ind,1)+2*real(obj.S_hat(z_ind,1)*exp(1i*x));
            plot_config.xlim_list=[1,0,max(obj.x_list)];
%             plot_config.xtick_list=[1,0,pi/2,pi,3*pi/2,2*pi];
%             plot_config.xticklabels_list={1,'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'};
            plot_config.ylim_list=[1,0,1];
            %plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
            if (max(obj.z_list)-min(obj.z_list))>1
                plot_config.ytick_list=0;
            else
                plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
            end

            plot_config.yticklabels_list={0};
            plot_config.label_list={1,'$x$','$z$'};
            plot_config.contour_line=0;
            plot_config.colorbar=1;
            %             plot_config.colormap='bluewhitered';
            if strcmp(p.my.z_basis_mode,'Fourier')
                plot_config.zlim_list=0;
                plot_config.ztick_list=0;
            else
                plot_config.zlim_list=[1,0,1];
                plot_config.ztick_list=[1,0,0.2,0.4,0.6,0.8,1];
            end
            %plot_config.user_color_style_marker_list={'k-','b--'};
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','isocontour_S.png'];
            plot_config.print_size=[1,1000,900];
            plot_contour(data,plot_config);
            if obj.no_ylabel
                plot_config.label_list{3}='';
                plot_config.name=[obj.h5_name(1:end-3),'_HB_','isocontour_S_no_ylabel.png'];
                plot_contour(data,plot_config);
            end

            x=obj.x_list;
            z_ind=1:length(obj.z_list);
            y=obj.z_list(z_ind);
            data{2}.x=x; data{2}.y=y;
            %data{2}.y=obj.z_list;

            %Update 2022/06/02, add the streamfunction corresponding to large-scale
            %shear
            data{2}.z=obj.psi_2D;
            [data{1}.x,data{1}.y]=meshgrid(x,y);
            data{1}.z=NaN*ones(size(data{1}.x));
            plot_config.xlim_list=[1,0,max(obj.x_list)];
%             plot_config.xtick_list=[1,0,pi/2,pi,3*pi/2,2*pi];
%             plot_config.xticklabels_list={1,'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'};
            plot_config.ylim_list=[1,0,1];
            plot_config.label_list={1,'$x$','$z$'};
            plot_config.streamline=2;
            plot_config.user_color_style_marker_list={'k-','r--'};
            plot_config.panel_num=2;
            plot_config.arrow_ratio=0.8;
            plot_config.linewidth=3;
            plot_config.colorbar=0;
%             plot_config.visible=1;
            plot_config.fontsize=28;
            plot_config.print_size=[1,500,900];
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','streamline.png'];
            plot_contour(data,plot_config);
            if obj.no_ylabel
                plot_config.label_list{3}='';
                plot_config.name=[obj.h5_name(1:end-3),'_HB_','streamline_no_ylabel.png'];
                plot_contour(data,plot_config);
            end

            
            %also output the obj to p.my.....
            p.my.obj=obj;
            p_old.my.folder_name=[p_old.my.folder_name,'/',plot_config.branch_name,'/graph_',...
                plot_config.point_name,'/'];
            p.my.folder_name=p_old.my.folder_name;
            disp_par(p);
        end
            
            
    case '2D_u_w_T_S_p'
        u_vel=p.u(1:p.np);
        w_vel=p.u(1+p.np:2*p.np);
        figure(101);
        u_vel_2D=reshape(u_vel,p.nz,p.nx);
        mesh(u_vel_2D);
        figure(102);
        plot(w_vel);
        figure(103);
        plot(mean(u_vel_2D,2));
    otherwise
            error('Wrong p.my.variable_version');
end





