% 2023-04-21 15:01:06.722925171 +0200
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


c = [1,1,0;
     0,0.5,0];

c = [1,1,0.25;
     0.25,0.5,0.25];


rng(1)

% isotropic pattern
b=generate_isotropic_pattern(4,200,1,0,0,0.1);
z=0.5*double(b>0.6);

figure(1)
clf
 imagesc(z);
 hold on;
 axis equal;
 axis off;
 axis tight;
 %colormap([1,1,1; 0.25,0.5,0.25]);
 colormap(c)
% colormap(flipud(gray));
if (pflag)
	pdfprint(1,'img2/pattern-spotted-schematic.pdf')
end
%x = linspace(0,1,50); z = cos(2*pi*5*x+0.*x'); imagesc(0.5*double(z>=0)); colormap(flipud(gray)); axis off; axis square; caxis([0,1]); pdfprint(1,'img/pattern-striped-schematic.pdf')
%x = linspace(0,1,50); z = cos(2*pi*5*x+0.*x'); imagesc(0.5*double(z>=0)); colormap([1,1,1;0,0.5,0]); axis off; axis square; caxis([0,1]); pdfprint(1,'img/pattern-striped-schematic.pdf')       
%x = linspace(0,1,50); z = cos(2*pi*5*x+0.*x'); imagesc(0.5*double(z>=0)); colormap(flipud(gray)); axis off; axis square; caxis([0,1]); pdfprint(1,'img/pattern-striped-schematic.pdf')         

if (1)
rng(100);
fc = 5;
lc = 1/fc;
Sc = 2.5;
L = 1;
n = 200;
x = linspace(0,1,n);
y = x';
[fx,fy,fr] = fourier_axis_2d(L*[1,1],n*[1,1]);
[a,b] = gampdf_mode2par(fc,Sc*lc);
Sx = gampdf(abs(fx),a,b);
%[a,b] = logn_mode2par(fc,Sc*lc);
%Sx = lognpdf(abs(fx),a,b);
%[a,b] = gamma_mode2par(1e-7*fc,Sc/lc);
Sy = exppdf(abs(fy),0.2/lc/Sc);
%gampdf(abs(fy),a,b);
 
T = sqrt(cvec(Sy)*rvec(Sx));

e = randn(n);

b = real(ifft2(sqrt(T).*fft2(e)));
if (1)
figure(20)
clf
subplot(2,2,1)
plot(fx/fc,Sx*fc)
hold on
plot(fy/fc,Sy*fc)
subplot(2,2,2)
end
figure(2)
clf
imagesc(b>quantile(b,0.6,'all'))
axis equal
axis off
 colormap(c)
if (pflag)
pdfprint(2,'img2/pattern-striped-schematic.pdf')
end
end

