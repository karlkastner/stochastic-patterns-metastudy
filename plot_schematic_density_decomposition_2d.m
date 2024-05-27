% Tue 11 Apr 11:51:28 CEST 2023
if (~exist('pflag','var'))
pflag = 0;
end
fflag = pflag;
if (1)
n=100;
%clf
%subplot(2,2,1)
splitfigure([2,2],[1,1],fflag);
cla
 fx = fftshift(fourier_axis(n/4.5,n));
 [a,b] = normpdf_mode2par(0,1);
 Sy=normpdf(fx,a,b);
 Sy=Sy/max(Sy);
 [a,b] = logn_mode2par(1,1);
 Sx=lognpdf(abs(fx),a,b);
 Sx=Sx/max(Sx);
 S=Sx*Sy';
 surface((fx),(fx),0.*(S),S,'edgecolor','none');
% view(-45,60);

 hold on;
fdx = fx>=0;
 plot3(min(fx)*ones(n/2,1),fx(fdx),Sx(fdx),'b','linewidth',1);
 plot3(min(fx)*ones(n/2,1),fx(~fdx),Sx(~fdx),':b','linewidth',1);
 plot3(fx(fdx),min(fx)*ones(n/2,1),Sy(fdx),'r','linewidth',1);
 plot3(fx(~fdx),min(fx)*ones(n/2,1),Sy(~fdx),':r','linewidth',1);
 grid on;
 xlabel('k_x/k_c');
 ylabel('k_y/k_c');
 zlabel('S/\lambda_c');
 colormap(flipud(gray));
 xlim(limits(fx));
 ylim(limits(fx)) 
 set(gca,'xdir','reverse')
 set(gca,'ydir','reverse')
av=5*4.5/6;daspect([av,av,1]); view([-45,30]);
shading interp
set(gcf,'Renderer','Painter')
gca_ = gca;
gca_.XLabel.Rotation = 30;
gca_.YLabel.Rotation = -30;
% axis square
end
if (1)
 fx = fftshift(fourier_axis(n/4.5,n)*n/(n-2));
 fdx = (fx>=0);
 k = mises_max2par(1.5);
 fr = hypot(fx,fx');
 Sr = lognpdf(fr,a,b);
 t=atan2(fx,fx');
 St = misesnpdf(t,0,k,6);
 S=Sr.*St;
%subplot(2,2,4)
splitfigure([2,2],[1,3],fflag);
cla
 surface(fx,fx,0.*S,S,'edgecolor','none');
% view(-45,60);
 xlabel('k_x/k_c');
 ylabel('k_y/k_c');
 zlabel('S/\lambda_c');
 av=15;daspect([av,av,1]); view([-45,30]);
 set(gca,'xdir','reverse')
 set(gca,'ydir','reverse')
shading interp

 colormap(flipud(gray));
axis tight;
 fr = fx(fdx);
 Sr = lognpdf(fr,a,b);
 t=linspace(-pi,pi,n)';
 St=misesnpdf(t,0,k,6);
%axis square
%subplot(2,2,4)
set(gcf,'Renderer','Painter')
gca_ = gca;
gca_.XLabel.Rotation = 30;
gca_.YLabel.Rotation = -30;
zlim([0,1])
%av=2.75;
%daspect([av,av*2/2,1]); view([-45,30]);
av=5*4.5/6;daspect([av,av,1]); view([-45,30]);
grid on
splitfigure([2,2],[1,4],fflag);
cla
 surface(t,fr,0.*(Sr*St'),4*Sr*St'.*(t'>=-pi/2&t'<=pi/2),'edgecolor','none');
 %view(-45,60);
 av=2.75;daspect([av,av*2/pi,1]); view([-45,30]);
 axis tight;
 hold on;
 plot3(pi/2*ones(n/2,1),fr,Sr,'b','linewidth',1);
 plot3(t,max(fr)*ones(n,1),4*St,'r','linewidth',1);
 ylim([0,max(fr)]);
%axis square
shading interp
xlim([-pi/2,pi/2])
set(gca,'xtick',pi*(-1/2:0.25:1/2),'xticklabel',{'-\pi/2','-\pi/4','0','\pi/4','\pi/2'},'xticklabelrot',0);
xlabel('\theta'); ylabel('k_r/k_c'); zlabel('S/\lambda_c'); grid on;
 colormap(flipud(gray));
set(gcf,'Renderer','Painter')
gca_ = gca;
gca_.XLabel.Rotation = 30;
gca_.YLabel.Rotation = -30;
end
if (pflag)
ps = 3.5;
pdfprint(11,'img/3d-anisotropic-decomposition.pdf',ps);
pdfprint(13,'img/3d-isotropic-x-y.pdf',ps);
pdfprint(14,'img/3d-isotropic-r-t.pdf',ps);
end
