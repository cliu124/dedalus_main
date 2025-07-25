clear all;
close all;
global mu alpha F RS1 tau;

RS1=0.01;
tau=0.01;
Pr=7;
alpha=1;

F=1/alpha^2*(1/Pr*alpha^4-1+(1+alpha^4)/tau);

xmesh=linspace(0,1,400);
solinit=bvpinit(xmesh,@guess);
sol=bvp4c(@bvpfcn_DDC,@bcfcn_DDC,solinit);


function dz_dtheta = bvpfcn_DDC(z,theta)
global mu alpha F RS1 tau;
dz_dtheta = zeros(size(theta));
dz_dtheta = [theta(2)
            1/3/alpha^2*(-RS1/tau*theta(1)-theta(1)*theta(4)+0*(1+alpha^4)*theta(1)*theta(6))
            theta(4)
            1/alpha^2*theta(1)*theta(2)
            theta(5)
            1/alpha^2/tau*theta(1)*theta(2)];
end

function res = bcfcn_DDC(ya,yb)
res = [ya(1)
       yb(1)
       ya(3)
       yb(3)
       ya(5)
       yb(5)];
end


function g = guess(x_num)
global mu alpha F RS1;

g = [sin(pi*x_num)
     cos(pi*x_num)
     sin(2*pi*x_num)
     cos(2*pi*x_num)
     sin(2*pi*x_num)
     cos(2*pi*x_num)
     ];
end
