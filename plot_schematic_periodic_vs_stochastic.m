% Mon  7 Aug 12:22:53 CEST 2023
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
qq = 1;
p=2;
 q=1/3;
 a=-log(1/4);
 a=-log(1/2);
% a=-log(1/3);
 % 1-b = -a
 b =1+a;
 b = 3;
b = 1;
%a = 0;
 L=1000;
% L = 10.5;
 dx=1/100;
 n=L/dx+1;
 x = cvec(innerspace(-L/2,L/2,n));
 y = [	(1+abs(x).^qq)./(1+b*abs(x).^qq).*cos(2*pi*x), ...
	exp(-a*abs(x).^p).*cos(2*pi*x), ...
	exp(-pi*a*abs(x).^p), ...
	zeros(n,1) ...
	... (q+(1-q)*exp(-a/(1-q)*abs(x).^p/(p))) ...
	];
 y((n+1)/2,end) = 1;
splitfigure([2,2],[1,1],fflag);
cla
 plot(x,y,'linewidth',1);
 xlim([-sqrt(eps),5.25]);
 ylim([-1,1]);
 hold on
 %hline(0) 
 area([3.25,4.25],[1,1]*0.99,-20,'facecolor','white','edgecolor','none');
 S=real(ifft(ifftshift(y,1)));
 df=1/L;
 S=2*S./(df*sum(S));
 fx=fourier_axis(x);
 xlabel('Lag distance $x/\lambda_c$','interpreter','latex');
 ylabel('Autocorrelation $R_x$','interpreter','latex')
 %legend('Stochastic','Periodic')
 legend('Periodic','Regular','Irregular','White Noise')
 % text(0,5,texlabel('$\frac{1}{2} \frac{L}{\lambda_c}$'));
 text(5,-1.1,'$\displaystyle\frac{1}{2} \frac{L}{\lambda_c}$','interpreter','latex');  
%,'intereter','latex');
 cm = colormap_krb();
 cm(end+1,:) = [0.5,0,0.5];
 set(gca,'colororder',cm);
% colormap(cm([3,2,1],:));

splitfigure([2,2],[1,2],fflag);
cla
 plot(fftshift(fx),fftshift(S,1),'linewidth',1);
 hold on
 xlim([0,2]);
 ylim([-sqrt(eps),6]);
 area([0.01,1.99],[5,5],3.5,'facecolor','white','edgecolor','none');
 xlabel('Wavenumber $k_x/k_c$','interpreter','latex');
 ylabel('Density $S_x/\lambda_c$','interpreter','latex')
 text(0,5,'$0.3 \frac{L}{\lambda_c}$','interpreter','latex');  
 %colormap(cm([3,2,1],:));
 set(gca,'colororder',cm);

 set(gca,'ytick',0:10);
if (pflag)
pdfprint(11,'img/stoch-vs-periodic-schematic-autocorrelation.pdf',3.5);
pdfprint(12,'img/stoch-vs-periodic-schematic-density.pdf',3.5);
end
