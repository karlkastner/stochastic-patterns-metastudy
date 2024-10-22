% Tue 11 Apr 11:51:28 CEST 2023
% Karl Kästner, Berlin
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
if (~exist('pflag','var'))
pflag = 0;
end
fflag = pflag;
if (1)
pdfx_str = 'normal';
% characteristic frequency
fc = 1;
% density maxima
Sxpc = 1;
Syc  = Sxpc;
% number of points
n  = 101;
% spectral axis
fx = fftshift(fourier_axis(n/4.5,n));
fy = fx;

splitfigure([2,2],[1,1],fflag);
cla
switch (pdfx_str)
case {'normal'}
	% parameter
	[a,b] = normalmirroredpdf_mode2par(fc,0.5*Sxpc);
	% density
	Sx = normalmirroredpdf(fx,a,b);
case {'lognormal'}
	[a,b] = lognpdf_mode2par(fc,Sxpc);
	Sx = lognmirroredpdf(fx,a,b);
end
% density along positive halve axis
Sxp = 2*Sx.*(fx>=0);
[a,b] = normpdf_mode2par(0,Syc);
Sy   = normpdf(fx,a,b);

% two-dimensional densityu
Sxy = cvec(Sx)*rvec(Sy);
surface(fx,fy,0.*Sxy,Sxy,'edgecolor','none');
% view(-45,60);

 hold on;
fdx = fx>=0;
 plot3(min(fx)*ones(n,1),fx,Sx,':b','linewidth',1);
 plot3(min(fx)*ones(n,1),fx,Sxp,'b','linewidth',1);
 plot3(fx,min(fx)*ones(n,1),Sy,'r','linewidth',1);
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
 %k = mises_max2par(1.5);
 fr = hypot(fx,fx');
 Sr = lognpdf(fr,a,b);
 t=atan2(fx,fx');
 k = misesn_max2par(0.5,6);
 St = misesnpdf(t,0,k,6);
 Srt = Sr.*St;
%subplot(2,2,4)
splitfigure([2,2],[1,3],fflag);
cla
 surface(fx,fx,0.*Srt,Srt,'edgecolor','none');
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
 [a,b] = lognpdf_mode2par(fc,1);
 Sr = lognpdf(fr,a,b);
 t=linspace(-pi,pi,n)';
 St  = misesnpdf(t,0,k,6);
 Stp = 2*St.*(t>-pi/2 & t <= +pi/2); 
%axis square
%subplot(2,2,4)
set(gcf,'Renderer','Painter')
gca_ = gca;
gca_.XLabel.Rotation = 30;
gca_.YLabel.Rotation = -30;
zlim([0,1])
%av=2.75;
%daspect([av,av*2/2,1]); view([-45,30]);
av=5*4.5/6;
daspect([av,av,1]);
view([-45,30]);
grid on
splitfigure([2,2],[1,4],fflag);
cla
surface(t,fr,0.*(Sr*St'),Sr*St'.*(t'>=-pi/2&t'<=pi/2),'edgecolor','none');
%view(-45,60);
av=2.75;daspect([av,av*2/pi,1]); view([-45,30]);
axis tight;
hold on;
plot3(pi/2*ones(length(fr),1),fr,Sr,'b','linewidth',1);
plot3(t,max(fr)*ones(n,1),Stp,'r','linewidth',1);
plot3(t,max(fr)*ones(n,1),St,'r:','linewidth',1);
ylim([0,max(fr)]);
%axis square
shading interp
xlim([-pi/2,pi/2])
set(gca,'xtick',pi*(-1/2:0.25:1/2),'xticklabel',{'-\pi/2','-\pi/4','0','\pi/4','\pi/2'},'xticklabelrot',0);
xlabel('\theta');
ylabel('k_r/k_c');
zlabel('S/\lambda_c');
grid on;
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
