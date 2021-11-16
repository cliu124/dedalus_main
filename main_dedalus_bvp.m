clear all;
close all;

slurm_num={'12516527'};

for slurm_ind=length(slurm_num)
    h5_name=['C:\Data\dedalus\dedalus_',...
            slurm_num{slurm_ind},...
            '\analysis\analysis_s1.h5'];
%     h5_name=['C:\Data\dedalus\dedalus_',...
%             slurm_num{slurm_ind},...
%             '\data.h5'];
    h5disp(h5_name);
    obj.z_list=h5read(h5_name,'/scales/z/1.0');
    obj.w_hat=h5read(h5_name,'/tasks/w_hat');
    obj.p_hat=h5read(h5_name,'/tasks/p_hat');
    obj.T_hat=h5read(h5_name,'/tasks/T_hat');
    obj.d_T_hat=h5read(h5_name,'/tasks/d_T_hat');
    obj.S_hat=h5read(h5_name,'/tasks/S_hat');
    obj.d_S_hat=h5read(h5_name,'/tasks/d_S_hat');
    obj.T_0=h5read(h5_name,'/tasks/T_0');
    obj.d_T_0=h5read(h5_name,'/tasks/d_T_0');
    obj.S_0=h5read(h5_name,'/tasks/S_0');
    obj.d_S_0=h5read(h5_name,'/tasks/d_S_0');

end
