% Thu 25 Jul 11:24:16 CEST 2024
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

% part 1 : demonstrate the effect of a spatial extents for one pattern

% characteristic frequency
fc = 0.5;
% regularity
regx = 1;
regy = 2;
% density maxima
Sxc = regx/fc;
Syc = regy/fc;

% spatial extent
Lx = 80;
Ly = 5;
% corresponding square extent
Lx_tiled = sqrt(Lx*Ly);
% number of tiles alon the longer axis
mx = Lx/Lx_tiled;
% spatial discretization
dx = 0.1;
dy = dx;
% spectral resolution
dfx = 1/Lx;
dfy = 1/Ly;
% number of points
nx = Lx/dx;
ny = Ly/dy;

% axes in frequency space
fx = fourier_axis(Lx,nx);
fy = fourier_axis(Ly,ny);

% density parameters
[fx0,sx] = normalwrappedpdf_mode2par(fc,Sxc)
[fy0,sy] = normpdf_mode2par(0,Syc);

% densities along axes
Sx = 0.5*normalwrappedpdf(fx,fx0,sx);
Sy = normpdf(fy,fy0,sy);
% two dimensional density
Sxy = cvec(Sx)*rvec(Sy);

% generate pattern 
rng(2);
e = randn(nx,ny);
b = real(ifft2(sqrt(Sxy).*fft2(e)));
% note that the density of the thresholded pattern is slightly different from Sxy,
% but sufficiently close for illustration here
b =b>quantile(b,0.65,'all');
sd = std(b,[],'all');
% periodogram for the entire domain
hatSxy = 1/sd^2*(Lx*Ly)/(nx^2*ny^2)*abs(fft2(b-mean(b,'all'))).^2;
barSx  = sum(hatSxy,2)*dfy;
barSx_smooth = ifftshift(meanfilt1(fftshift(barSx),mx));

% number of points per tile
nx_ = Lx_tiled/dx;
% frequency axis per tile
fx_tiled = fourier_axis(Lx_tiled,nx_);
% split into tiles 
hatS_avg = 0;
for idx=1:mx
	% tile of pattern
	b_     = b((1:Lx_tiled/dx)+(idx-1)*Lx_tiled/dx,:);
	% its standard deviation
	sd_    = std(b,[],'all');
	% its periodogram
	hatS_i = 1/sd_.^2*Lx_tiled*Ly/(nx_^2*ny^2)*abs(fft2(b_-mean(b_,'all'))).^2;
	% stack periodograms
	hatS_avg  = hatS_avg + hatS_i;
end
hatS_avg = hatS_avg/mx; 
barSx_tiled  = sum(hatS_avg,2)*dfy;

% plot pattern
splitfigure([2,2],[1,1],fflag);
imagesc(b)
axis equal
axis tight
colormap(colormap_vegetation2())
axis off
hold on;
for idx=1:4;
	plot([0,Ly,Ly,0,0]/dy,[0,0,Lx_tiled,Lx_tiled,0]/dx+(idx-1)*Lx_tiled/dx,'r','linewidth',1);
end

% plot rearranged tiles of pattern
splitfigure([2,2],[1,2],fflag);
b_ = [];
for idx=1:4
	b_ = [b_,b((1:Lx_tiled/dx)+(idx-1)*Lx_tiled/dx,:)];
	%reshape(b,[],Lx_tiled/dx);
end
imagesc(b_)
axis equal
axis tight
colormap(colormap_vegetation2())
axis off
hold on;
for idx=1:4;
 plot([0,Ly,Ly,0,0]/dy+(idx-1)*Ly/dy,[0,0,Lx_tiled,Lx_tiled,0]/dx+0*(idx-1)*Lx_tiled/dx,'r','linewidth',1);
end

% plot density
splitfigure([2,2],[1,3],fflag);
cla
fdx = fx>=0;
plot(fx(fdx)/fc,2*fc*barSx(fdx),'.');
hold on
plot(fx_tiled/fc,2*barSx_tiled*fc,'r.')
plot(fx(fdx)/fc,2*fc*barSx_smooth(fdx),'r','linewidth',1)
plot(fx(fdx)/fc,2*fc*Sx(fdx),'b','linewidth',1)
xlim([0,3])
%set(gca,'colororder',(colormap_krb))
xlabel('Wavenumber k/k_c')
ylabel('Density S_x/\lambda_c')
axis square
legend('Periodogram','Tiled','Smoothed','Density')

if (pflag)
	ps = 3.5;
	pdfprint(11,'img/pattern-rectangular.pdf',ps);
	pdfprint(12,'img/pattern-rectangular-rearranged.pdf',ps);
	pdfprint(13,'img/pattern-rectangular-density.pdf',ps);
end

% part 2 : estimate standard error for spatial extents with different aspect ratio
if (1)

rng(0)
nrep = 100;
Lxmax = 100;
Lx  = Lxmax*(logspace(-1,0,10));
%Ly = 100./Lx;
regx = 1;
regy = 1;
fc   = 1;
Sxc  = regx/fc;
Syc  = regy/fc;
[fx0,sx] = normalwrappedpdf_mode2par(fc,Sxc)
[fy0,sy] = normpdf_mode2par(0,Syc);
dx   = 0.1/fc;

regx_est = [];
for idx=1:length(Lx)
	idx
	Ly = Lxmax/Lx(idx);
	nx = round(Lx(idx)/dx);
	ny = round(Ly/dx);
	mx = sqrt(Lx(idx)/Ly);
	fx = fourier_axis(Lx(idx),nx);
	fy = fourier_axis(Ly,ny);
	Sx  = 0.5*normalwrappedpdf(fx,fx0,sx);
	Sy = normpdf(fy,fy0,sy);
	Sxy = cvec(Sx)*rvec(Sy);
	Txy = sqrt(Sxy); 
for jdx=1:(nrep)
	e = randn(nx,ny);
	b = ifft2(Txy.*fft2(e));
	b = real(b); 
	% periodogram
	hatSxy = abs(fft2(b-mean(b,'all'))).^2;
	% normalize
	hatSxy = hatSxy/(sum(hatSxy)/Lx(idx)/Ly);
	% density component
	barSx = sum(hatSxy,2)/Ly;
	% smoothed density component
	mx_(idx) = round(mx);
	barSx_ = ifftshift(meanfilt1(fftshift(barSx),mx_(idx)));

	[Sxc,mdx] = max(barSx);
	fc_        = fx(mdx);
   	regx_est(idx,jdx,1) = 2*Sxc*fc_;

	[Sxc,mdx] = max(barSx_);
	fc_        = fx(mdx);
   	regx_est(idx,jdx,2) = 2*Sxc*fc_;
end % for jdx
end % for idx
res  = regx_est-regx;
bias = squeeze(mean(res,2));
sd   = squeeze(std(res,[],2));
rmse = squeeze(rms(res,2));

Ly = Lxmax./Lx;
splitfigure([2,3],[2,1],fflag);
semilogx(Lx./Ly,rmse/regx);
xlabel('Aspect ratio $L_x/L_y$','interpreter','latex');
ylabel({'Relative root mean square','error $\displaystyle \frac{||\overline \mathrm{reg}_x - \mathrm{reg}_x||_2}{\mathrm{reg}_x}$'},'interpreter','latex');
legend('Unsmoothed','Smoothed','location','northwest');
axis square

splitfigure([2,3],[2,2],fflag);
semilogx(Lx./Ly,bias);
ylabel('Bias');

splitfigure([2,3],[2,3],fflag);
semilogx(Lx./Ly,sd);
ylabel('sd');

if (pflag)
	ps = 3.5;
	pdfprint(21,'img/pattern-rectangular-error.pdf',ps);
end

end
