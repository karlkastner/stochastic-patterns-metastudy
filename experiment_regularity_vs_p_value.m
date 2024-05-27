% Fri 17 Feb 14:45:58 CET 2023
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
%% numerical experiment, demonstrating the unsuitability of p-values as gradual
%% regularity estimates
%
if (~exist('pflag','var'))
	pflag = 0;
end
fflag = pflag;

if (~exist('qr','var'))
% reset random nummer generator for exact reproducibility
rng(0);
% number of samples to estimate test quantiles  
np = 100;
% spatial extend of domain
L  = 10;
% spectral resolution
df    = 1./L;
% central frequency fc = 1/lambda_c
fc = 1;
% spatial resolution
dx = 1/(fc*50);
% domain size in pixels
n  = round(L/dx);
% smoothing radius for estimating the density in bins
nf = sqrt(L*fc);
% axes in frequency space
[fx,fy] = fourier_axis_2d([L,L],[n,n]);
% number of samples of monte-carlo estimate of test
% not necessary to set, as spatial domain is not masked
ns = 101;
% frequency range to include into periodicity test
% measured in bins containing the upper pfmks-fraction of spectral energy
pfmsk = 0.9;

% repeat experiment for a range of regularities
regularity = 2.^(-3:3);

% quantiles of estimates
qp = zeros(3,length(regularity));
qr = zeros(3,length(regularity));
for jdx=1:length(regularity)
	disp(jdx);
	% max of density
	Sc      = regularity(jdx)/fc;
	% generate log-normal density for Sx
	[a,b]   = logn_mode2param(fc,Sc);
	Sx      = lognpdf(abs(fx),a,b); 
	% generate exp-density for Sy
	c       = exppdf_max2par(Sc);
	Sy      = exppdf(abs(fy),c);
	% normalize
	Sx    = 2*Sx/(sum(Sx)*df);
	Sy    = 2*Sy/(sum(Sy)*df);
	% transfer function	
	T     = sqrt(0.5*Sx)*sqrt(0.5*Sy');

	% sample np repetitions
	hat_reg = zeros(np,1);
	p = zeros(np,1);
	for idx=1:np
		% white noise
		e    = randn(n);
		% patterns
		b    = ifft2(T.*fft2(e));
		% periodogram
		hatS = abs(fft2(b-mean(b(:)))).^2;
		barS = gaussfilt2(hatS,nf);
		% mask frequency space, points containing 90% of spectral energy
		[sS,sds] = sort(barS(:),'descend');
		iS   = cumsum(sS(:));
		iS   = iS/iS(end);
		fdx  = find(iS>pfmsk,1,'first');
		fmsk = false(n);
		fmsk(sds(1:fdx)) = true;
		% exclude symmetric part from text
		fmsk(fx<0) = 0;
		bmsk = [];
		% test for periodicity
		[issignificant,p(idx)] = periodogram_test_periodicity_2d(b,nf,bmsk,fmsk,ns);
		
		% estimate density along x
		hatSx  = sum(hatS,2)*df;
		% normalize
		hatSx  = 2*hatSx/(sum(hatSx)*df);
		% extract regularity
		[hat_Sc,mdx] = max(hatSx);
		hat_fc = abs(fx(mdx));
		hat_reg(idx) = hat_Sc.*hat_fc;
	end % for jdx
	% quartiles
	qp(:,jdx) = quantile(p,[0.25,0.5,0.75]);
	qr(:,jdx) = quantile(hat_reg,[0.25,0.5,0.75]);
end % idx
end % if not exist qr

% plot p-values
lim = [0.7,1.4].*limits(regularity)';
splitfigure([2,2],[1,1],fflag);
cla();
errorbar(regularity,qp(2,:),qp(2,:)-qp(1,:),qp(3,:)-qp(2,:),'*');
xlim(lim);
ylim([0,1])
set(gca,'xscale','log');
xlabel('True Regularity S_{cx}/\lambda_c');
ylabel('p');
set(gca,'xtick',2.^(-3:3))

% plot regularity
splitfigure([2,2],[1,2],fflag);
cla();
errorbar(regularity,qr(2,:),qr(2,:)-qr(1,:),qr(3,:)-qr(2,:),'*'); xlim(lim); set(gca,'xscale','log','yscale','log');
hold on;
plot(lim,lim,'k-');
ylim(lim);
xlabel('True Regularity S_{cx}/\lambda_c');
ylabel('Estimated regularity S_{cx}/\lambda_c');
set(gca,'xtick',2.^(-3:3))
set(gca,'ytick',2.^(-3:3))

if (pflag)
	ps = 3.5;
	pdfprint(11,'img/p-vs-regularity.pdf',ps);
	pdfprint(12,'img/regularity-estimate-vs-regularity.pdf',ps);
end

