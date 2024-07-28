% Mon  1 Jul 12:14:00 CEST 2024
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
% estimate the error of the regularity estimate for various regularity measures
%
if (~exist('pflag','var'))
	pflag =0;
end
fflag = pflag;

% characteristic frequency
fc = 1;
% density maximum
Sc = 1;
% spatial extent
Lx  = 10/fc;
% spatial resolution
dx = 0.1*fc;

% number of points
nx = Lx/dx;

% spatial axis
x = (0:nx-1)'*Lx/nx;
y = x;
% spectral axis
fx = fourier_axis(x);
fy = fourier_axis(y);
% spectral resolution
dfx = 1./Lx;
dfy = 1./Lx;

% numer of repetition of monte-carlo analysis
m = 1e3;

% added noise levels
% turns out that the added noise  does not has a large influence
pe = [0,0.2];

distribution_C = {'gauss', 'laplace', 'cauchy'};
% for each distribution
for ddx=1:length(distribution_C)
	distribution = distribution_C{ddx};
	% for each added noise level
	for pdx=1:length(pe)
		
		rng(0)
		reg = [];
		reg_ = [];
		switch (distribution)
		case {'gauss'}
			% density parameter
			[f0x,sx] = normalwrappedpdf_mode2par(fc,Sc);
			[f0y,sy] = normpdf_mode2par(0,Sc);
			% density
			Sx = normalwrappedpdf(fx,f0x,sx);
			Sy = normpdf(fy,f0y,sy);
		case {'laplace'}
			[f0x,sx] = laplacewrappedpdf_mode2par(fc,Sc);
			Sx = laplacewrappedpdf(fx,f0x,sx);
			[f0y,sy] = laplacepdf_mode2par(0,Sc);
			Sy = laplacepdf(fy,f0y,sy);
		case {'cauchy'}
			[f0x,sx] = cauchywrappedpdf_mode2par(fc,Sc);
			Sx = cauchywrappedpdf(fx,f0x,sx);
			[f0y,sy] = cauchypdf_mode2par(0,Sc);
			Sy = cauchypdf(fy,f0y,sy);
		end
		% two dimensional density
		Sxy = cvec(Sx)*rvec(Sy);

		% estimate regularity from density
		reg_(1,:) = regularity_measure(fx,Sx,distribution);

		% monte-carlo iteration
		for idx=1:m
			% generate a random pattern
			e  = randn(nx);
			b = ifft2(sqrt(Sxy).*fft2(e));
			b = real(b);
			% standard deviation
			sd = std(b,[],'all');
			% add white noise
			ew = randn(nx);
			b = sqrt(1-pe(pdx)^2)*b + pe(pdx)*sd*ew;
		
			% periodogram
			hatS = abs(fft2(b-mean(b,'all'))).^2;
		
			% density estimates
			barSx = sum(hatS,2)*dfy;
			barSxp = barSx.*(fx>=0);
			% normalize
			barSxp = barSxp./(sum(barSxp)*dfx);
			
			% estimate regularity from patterns
			[reg(idx,:),leg_C] = regularity_measure(fx,barSx,distribution);
		end
		m = size(reg,2);
		
		splitfigure([2,3],[1,ddx + 3*(pdx-1)],fflag);
		cla
		q = quantile(reg,[0.16,0.5,0.84]);
		errorbar(1:m,q(2,:),q(1,:)-q(2,:),q(3,:)-q(2,:),'*');
		set(gca,'xtick',1:m,'xticklabel',leg_C);
		xlim([0.5,m+0.5]);
		ylim([0,3.5]);
		ylabel('Estimated regularity $S_x/\lambda_c$','interpreter','latex');
		hold on
		plot(1:m,reg_,'*r')
		title([upper(distribution(1)),distribution(2:end)]);
	end % for pdf
	sd_reg(ddx,:) = std(reg);
end % for ddx each distribution
sd_reg

if (pflag)
	ps = 3;
	aspect=1.2;
	pdfprint(11,'img/regularity-measure-uncertainty-gauss.pdf',ps,aspect);
	pdfprint(12,'img/regularity-measure-uncertainty-laplace.pdf',ps,aspect);
	pdfprint(13,'img/regularity-measure-uncertainty-cauchy.pdf',ps,aspect);
end

