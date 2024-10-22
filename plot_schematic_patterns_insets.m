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

% colormap
c = [1,1,0;
     0,0.5,0];

c = [1,1,0.25;
     0.25,0.5,0.25];

% reset random number generator for reproducibility
rng(1)

% isotropic pattern
%[b,x,y]=generate_isotropic_pattern(1/lc,n,L,alpha,[],[],0,s(idx));
fc = 4;
nx = 200;
L = 1;
alpha = 0;
p = 1;
q = 1;
scale = [];
st = 0.1;
b = generate_isotropic_pattern(fc,nx,L,alpha,p,q,[],st);
% threshold
z = 0.5*double(b>0.6);

figure(1)
clf()
imagesc(z);
hold on;
axis equal;
axis off;
axis tight;
colormap(c)

if (pflag)
	pdfprint(1,'img2/pattern-spotted-schematic.pdf')
end

if (1)
rng(100);
% characteristic wavelength
lc = 1;
% characteristic frequency
fc = 1/lc;
% regularity reg = Sc/lc = Sc fc
reg = 4;
% density maximum
Sxpc = reg/fc;
Syc  = Sxpc;
% spatial extent
Lx = 5*lc;
% spatial resolution
dx = lc/40;
% number of points
nx = Lx/dx;
%x = linspace(0,1,n);
%y = x';
[fx,fy,fr] = fourier_axis_2d(Lx*[1,1],nx*[1,1]);
[a,b] = gampdf_mode2par(fc,Sxpc/lc);
Sx    = 0.5*gampdf(abs(fx),a,b);
%[a,b] = logn_mode2par(fc,Sc/lc);
%Sx = lognpdf(abs(fx),a,b);
%[a,b] = gamma_mode2par(1e-7*fc,Sc/lc);
[fy0,sy] = laplacepdf_mode2par(0,Syc)
Sy = laplacepdf(fy,fy0,sy);
%Sy = exppdf(abs(fy),0.2*Sxpc);
%gampdf(abs(fy),a,b);

% transfer function 
T = sqrt(cvec(Sy)*rvec(Sx));
% generate random pattern
e = randn(nx);
b = real(ifft2(sqrt(T).*fft2(e)));
figure(20)
clf
subplot(2,2,1)
plot(fx/fc,Sx*fc)
hold on
plot(fy/fc,Sy*fc)
subplot(2,2,2)

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

