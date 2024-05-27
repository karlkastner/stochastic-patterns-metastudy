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
%% numerical experiment demonstrating the bias and standard error of the
%% estimated regularity
%

if (~exist('pflag','var'))
	pflag = 0;
end
fflag = pflag;

% recompute the experiment only if it has not yet run
if (~exist('serr_rel','var'))
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
dx = 1/(50*fc);

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
	% construct density
	n     = round(L(jdx)/dx);
	Sc    = regularity(idx)/fc;
	[a,b] = logn_mode2param(fc,Sc);
	c     = exppdf_max2par(Sc);
	[fx,fy] = fourier_axis_2d([L(jdx),L(jdx)],[n,n]);
	% generate log-normal and exp density
	Sx    = lognpdf(abs(fx),a,b); 
	Sy    = exppdf(abs(fy),c);
	% spectral resolution
	df = 1./L(jdx);
	% normalize
	Sx    = 2*Sx/(sum(Sx)*df);
	Sy    = 2*Sy/(sum(Sy)*df);
	% transfer function
	T     = sqrt(0.5*Sx)*sqrt(0.5*Sy');
	% repeat exeriment for estimating the bias and variation
	hat_Sc = zeros(m,1);
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
		hatSx = 2*hatSx/(sum(hatSx)*df);
		% estimate regularity
		[hat_Sc(kdx),mdx] = max(hatSx);
		hat_fc(kdx)       = abs(fx(mdx));
	end % for kdx
	hat_reg = hat_Sc.*hat_fc;
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
xlabel('Regularity $S_{cx}/\lambda_c$','interpreter','latex');
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
xlabel('Regularity $S_{cx}/\lambda_c$','interpreter','latex');
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
if (~pflag)
	title('Relative Bias');
end
xlabel('Regularity $S_{cx}/\lambda_c$','interpreter','latex');
ylabel('Spatial extent $L/\lambda_c$','interpreter','latex');

if (pflag)
	ps = 2;
	pdfprint(11,'img/regularity-estimate-serr.pdf',ps);
	pdfprint(12,'img/regularity-estimate-sd.pdf',ps);
	pdfprint(13,'img/regularity-estimate-bias.pdf',ps);
end

