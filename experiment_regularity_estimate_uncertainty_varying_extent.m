% Fri 17 Feb 12:48:52 CET 2023
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
%
%% numerical experiment estimating the bias and standard error of the
%% regularity estimated based on the height of the mode
%% spatial extend is varied
%

if (~exist('pflag','var'))
	pflag = 0;
end
fflag = pflag;

% recompute the experiment only if it has not yet run
if (~exist('serr_rel','var'))
pdfx_str = 'lognormal';
pdfy_str = 'normal';

% reset random number generator for exact reproducibility
rng(0)
% number of samples
m  = 100;
%regularity = logspace(-1,1,5);
regularity = 2.^(-3:4);
% spatial extend of domain
L          = 2.^(1.5:0.5:4.5);
% characteristic frequency fc = 1/lambda_c
fc         = 1;
% spatial resolution
dx = 1/(20*fc);

nf = 0;

% allocate space
bias     = zeros(length(regularity),length(L));
bias_rel = zeros(length(regularity),length(L));
serr     = zeros(length(regularity),length(L));
serr_rel = zeros(length(regularity),length(L));
sd       = zeros(length(regularity),length(L));
sd_rel   = zeros(length(regularity),length(L));
for idx=1:length(regularity)
% display progress
disp(idx)
for jdx=1:length(L)
	n       = round(L(jdx)/dx);
	[fx,fy] = fourier_axis_2d([L(jdx),L(jdx)],[n,n]);
	% construct density
	Sxpc    = regularity(idx)/fc;
	Syc     = Sxpc;
	[a,b]   = lognmirroredpdf_mode2par(fc,0.5*Sxpc);
	% generate mirrored log-normal density
	Sx    = lognmirroredpdf(fx,a,b);

	switch (pdfy_str)
	case {'laplace'}
		% note that for laplace, the outer product Sx Sy ceases to have elliptic contours
		c     = laplacepdf_max2par(Syc);
		Sy    = laplacepdf(fy,0,c);
	case {'normal'}
		[f0,s] = normpdf_mode2par(0,Syc);
		Sy     = normpdf(f0,s);
	end
	% spectral resolution
	df = 1./L(jdx);
	% normalize
	Sx    = Sx/(sum(Sx)*df);
	Sy    = Sy/(sum(Sy)*df);
	% transfer function
	T     = sqrt(Sx)*sqrt(Sy');
	% repeat exeriment for estimating the bias and variation
	hat_Sxpc = zeros(m,1);
	hat_fc = zeros(m,1);
	for kdx=1:m
		% white noise
		e  = randn(n);
		% pattern
		b  = ifft2(T.*fft2(e));
		% periodogram (not normalized)
		hatS  = abs(fft2(b-mean(b,'all'))).^2;
		% estimate density
		hatSx = sum(hatS,2)*df;
		hatSxp = hatSx.*(fx>=0);
		hatSxp = hatSxp/(sum(hatSxp)*df);
		% estimate regularity
		[hat_Sxpc(kdx),mdx] = max(hatSxp);
		hat_fc(kdx)       = abs(fx(mdx));
	end % for kdx
	hat_reg = hat_Sxpc.*hat_fc;
	% bias
	bias(idx,jdx) = mean(hat_reg)-regularity(idx);
	bias_rel(idx,jdx) = bias(idx,jdx)./regularity(idx);
	% standard deviation
	sd(idx,jdx)     = std(hat_reg);
	sd_rel(idx,jdx) = sd(idx,jdx)/regularity(idx);
	% standard error
	serr(idx,jdx) = rms(hat_reg - regularity(idx));
	serr_rel(idx,jdx) = serr(idx,jdx)./regularity(idx);
	end % for jdx
end % for idx

end % if not yet run

% plot standard error
splitfigure([2,3],[1,1],fflag);
cla();
contourf(regularity,L,serr_rel');
shading interp
xlabel('Regularity $S_{xc}^+/\lambda_c$','interpreter','latex');
ylabel('Spatial extent $L/\lambda_c$','interpreter','latex');
set(gca,'xscale','log','yscale','log')
axis square
colorbar('location','southoutside');
set(gca,'xtick',2.^(-3:5));
set(gca,'ytick',2.^(-3:5));
shading interp
if (~pflag)
	title('Relative Serr');
end

% plot standard deviation
splitfigure([2,3],[1,2],fflag);
cla();
contourf(regularity,L,sd_rel');
set(gca,'xscale','log','yscale','log')
colorbar('location','southoutside');
axis square
shading interp
set(gca,'xtick',2.^(-3:5));
set(gca,'ytick',2.^(-3:5));
if (~pflag)
	title('Relative Standard Deviation');
end
%xlabel('Regularity $S_{cx}/\lambda_c$','interpreter','latex');
xlabel('Regularity $S_{xc}^+/\lambda_c$','interpreter','latex');
ylabel('Spatial extent $L/\lambda_c$','interpreter','latex');

% plot bias
splitfigure([2,3],[1,3],fflag);
cla();
contourf(regularity,L,bias_rel');
set(gca,'xscale','log','yscale','log')
axis square
colorbar('location','southoutside');
shading interp
set(gca,'xtick',2.^(-3:5));
set(gca,'ytick',2.^(-3:5));
clim = max(abs(bias_rel(:)))*[-1,1];
caxis(clim);
if (~pflag)
	title('Relative Bias');
end
%xlabel('Regularity $S_{cx}/\lambda_c$','interpreter','latex');
xlabel('Regularity $S_{xc}^+/\lambda_c$','interpreter','latex');
ylabel('Spatial extent $L/\lambda_c$','interpreter','latex');

if (pflag)
	ps = 2;
	base = sprintf('img/regularity-estimate-Sx-%s-Sy-%s',pdfx_str,pdfy_str);
	pdfprint(11,[base,'-serr.pdf'],ps);
	pdfprint(12,[base,'-sd.pdf'],ps);
	pdfprint(13,[base,'-bias.pdf'],ps);
end

