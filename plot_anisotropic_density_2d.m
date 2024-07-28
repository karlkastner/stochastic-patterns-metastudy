% Mon 22 Jul 10:34:01 CEST 2024
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
% l=normcdf(-1.5:0.5:1.5);
 l=innerspace(0,1,11)
fx = linspace(0,2.5);
 fy = linspace(-1.25,1.25,101);
 Sxpc = 2.^(-0.5:0.5:0.5);
 Syc = Sxpc;
% p=2.^[-1/2,0,1/2];
 for idx=1:3;
 for jdx=1:3
 [a,b]=normalwrappedpdf_mode2par(1,Sxpc(idx)); %R*p(idx));
 Sx = 0.5*normalwrappedpdf(fx,a,b);
 [a,b]=normpdf_mode2par(0,Syc(jdx)); %c/p(idx));
 Sy=normpdf(fy,a,b);
% subplot(1,3,idx);
splitfigure([3,3],[1,(idx-1)*3+jdx],fflag);
 Sxy = (cvec(Sx)*rvec(Sy))';
 [c,h]=contourf(fx,fy,Sxy,l);
 axis equal;
if (0)
 xlabel('k_x/k_c');
 ylabel('k_y/k_c');
else
	set(gca,'xticklabel',{});
	set(gca,'yticklabel',{});
end
 colormap(flipud(colormap('gray')))
 %text(0.5,1,sprintf('$$S_{xc}^+ \\cdot S_{yc} = %1.2g \\frac{S_{yc}}{S_{xc}^+} = %1.2g$$',Sxc*Syc,1./p(idx).^2),'interpreter','latex');
 caxis([0,1]);
if (idx<3)
% text(0.05,0.95,sprintf('$$S_{xc}^+ \\cdot S_{yc} = %1.2g$$',Sxc*Syc),'interpreter','latex');
% text(1.545,0.92,sprintf('$$\\frac{S_{yc}}{S_{xc}^+} = %1.2g$$',1./p(idx).^2),'interpreter','latex');
else
% text(1.5,0.95,sprintf('$$S_{xc}^+ S_{yc} = %1.2g$$',Sxc*Syc),'interpreter','latex');
% text(1.545,0.92,sprintf('$$\\frac{S_{yc}}{S_{xc}^+}=\\frac{1}{2}$$',1./p(idx).^2),'interpreter','latex');
end
%axis off
if (pflag)
	ps = 4;
	name = sprintf('img/density-Sxy-Sxpc-%1.2g-Syc-%1.2g.pdf',Sxpc(idx),Syc(jdx));
	pdfprint(10+(idx-1)*3 + jdx,name,ps);
end
 end;
end

