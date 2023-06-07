clear all;
close all;
clc;

%Ra_{T,q}=11400, Pr=0.1, slurm_num={'15376325'};
load('quasi_periodic_11600.mat');

% ind_local_max=islocalmax(data{1}.y);
freq_peak=freq_local_max';
% f1=(freq_peak(1)+freq_peak(2))/2;
% f2=(freq_peak(2)-freq_peak(1))/2;%freq_peak(27)-freq_peak(1)*23;

f1=freq_peak(1);
f2=freq_peak(2);
% 
% f1=0.1482;
% f2=0.2267;
freq=[];
for m=0:10
    for n=0:10
        freq=[freq,m*f1+n*f2];
    end
end
freq=sort(freq)';

error('1');

plot_config.label_list={1,'$f $','PSD'};
            
plot_config.name='PSD_quasi_periodic.png';

%plot_config.xlim_list=[1,0,10];
plot_config.print=1;
plot_config.visible=1;
plot_line(data,plot_config);
