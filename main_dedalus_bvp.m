clear all;
close all;

slurm_num={'12514099'};

for slurm_ind=length(slurm_num)
    h5_name=['C:\Data\dedalus\dedalus_',...
            slurm_num{slurm_ind},...
            '\analysis\analysis_s1.h5'];
%     h5_name=['C:\Data\dedalus\dedalus_',...
%             slurm_num{slurm_ind},...
%             '\data.h5'];
    h5disp(h5_name);
    obj.z=h5read(h5_name,'/z');
    obj.w_hat=h5read(h5_name,'/w_hat');
    obj.p_hat=h5read(h5_name,'/p_hat');
    obj.T_hat=h5read(h5_name,'/T_hat');
    obj.d_T_hat=h5read(h5_name,'/d_T_hat');
    obj.S_hat=h5read(h5_name,'/S_hat');
    obj.d_S_hat=h5read(h5_name,'/d_S_hat');
    obj.T_0=h5read(h5_name,'/T_0');
    obj.d_T_0=h5read(h5_name,'/d_T_0');
    obj.S_0=h5read(h5_name,'/S_0');
    obj.d_S_0=h5read(h5_name,'/d_S_0');

end
