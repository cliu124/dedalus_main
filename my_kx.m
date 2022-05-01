function [kx, ky]=my_kx(p,u)

par=u(p.nu+1:end); 
kx=par(1);
ky=par(2);
Ra_T=par(3);
Ra_S2T=par(4);
Le=par(5);
Pr=par(6);
dy_T_mean=par(7);
dy_S_mean=par(8);

tau=1/Le;
Ra_S=Ra_S2T/tau;
R_rho_T2S=Ra_T/Ra_S2T;
if kx==-10^5
    %The option to use the wavenumber scaling from Yang.
    kx=2*pi/(2*14.8211*Ra_S^(-0.2428)/R_rho_T2S^(0.25/2));
elseif kx==-10^3
    kx=sign(ky)*sqrt(p.my.kx_square-ky^2);
end

if ky==-10^5
    %The option to use the scaling from Yang
    ky=2*pi/(2*14.8211*Ra_S^(-0.2428)/R_rho_T2S^(0.25/2));
elseif ky==-10^4
    %The option to set the ky the same as kx
    ky=kx;
elseif ky==-10^3
    ky=sign(kx)*sqrt(p.my.kx_square-kx^2);
end

end