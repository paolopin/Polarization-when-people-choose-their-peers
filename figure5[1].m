gg=grid_fine;
R=linspace(0.1,0.9,grid_fine);
for i=1:gg
    diam2(i)=(32*f(i)-36*f(i)^2+9*f(i)^3)/(9*f(i)^2)-(4/9)*((-(-64*f(i)^2+144*f(i)^3-108*f(i)^4+27*f(i)^5)/(f(i)^4)))^(1/2);
    diam3(i)=(1/36)*((f(i)*(-5+3*f(i))*(-4+3*f(i)))/(-1+f(i))^2-2*(-(f(i)^2*(-4+3*f(i))^3)/(-1+f(i))^4)^(1/2));
    diam4(i)=(f(i)*(-4+3*f(i))*(-40+27*f(i)))/(32-27*f(i))^2-12*(-(f(i)^2*(-4+3*f(i))^3)/(-32+27*f(i))^4)^(1/2);
    diam5(i)=(1/36)*(-4*(-((f(i)^2*(-4+3*f(i))^3/(-5+4*f(i))^4)))^(1/2)+(f(i)*(68+9*f(i)*(-11+4*f(i))))/(5-4*f(i))^2);
    diam6(i)=-(20/9)*(-(f(i)^2*(-4+3*f(i))^3)/(-32+25*f(i))^4)^(1/2)+((f(i)*(-4+3*f(i))*(-104+75*f(i)))/(96-75*f(i))^2);
    diam2a(i)=(4/9)*((-8*(-(f(i)^2*(-4+3*f(i))^3)/(4+f(i))^4)^(1/2))+((f(i)*(80+9*(-8+f(i))*f(i)))/(4+f(i))^2));
end
%figure
hold on
contourf(f,V,W')
colorbar

ylabel('V')
xlabel('f')

plot(f,diam2,'r','LineWidth',2)
plot(f,diam3,'y','LineWidth',3)
plot(f,diam4,'LineWidth',2)
plot(f,diam5,'LineWidth',2)
plot(f,diam6,'LineWidth',2)
hold off

