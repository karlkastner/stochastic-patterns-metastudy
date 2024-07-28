% Sun 30 Jun 13:06:49 CEST 2024
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
% demonstrate that the transformed regularity measures based on different density
% properties converge to the same value with increasing regularity,
% when the distribution, and hence transformation, is known
%
% estimates are conducted for the raw and the thresholded pattern

if (~exist('pflag','var'))
	pflag =0;
end
fflag = pflag;
rng(0);
% characteristic frequency
fc = 1;
% density maxima
Sxpc = logspace(-0.9,1,3e1);
% spatial extent
Lx = 100/fc;
% spatial resolution
dx = fc/100;
% number of points
nx = Lx/dx;
% spatial axis
x = (0:nx-1)'*Lx/nx;
% spectral axis
fx = fourier_axis(x);
% number of patterns for monte carlo analysis
m = 1000;
e = randn(nx,m);

distribution_C = {'gauss','laplace','cauchy'}

clear out out_thresh
for jdx=1:length(distribution_C)
distribution = distribution_C{jdx};

reg = [];
reg_thresh = [];
for idx=1:length(Sxpc)
switch (distribution)
case {'gauss'}
	% density parameter
	[f0,s] = normalwrappedpdf_mode2par(fc,Sxpc(idx));
	% density
	S  = normalwrappedpdf(fx,f0,s);
	% width
	w = normalwrappedpdf_width(f0,s);
case {'laplace'}
	[f0,s] = laplacewrappedpdf_mode2par(fc,Sxpc(idx));
	S = laplacewrappedpdf(fx,f0,s);
	w = laplacewrappedpdf_width(f0,s);
case {'cauchy'}
	[f0,s] = cauchywrappedpdf_mode2par(fc,Sxpc(idx));
	f0_(idx,1) = f0;
	s_(idx,1) = s;
	S = cauchywrappedpdf(fx,f0,s);
	try
	w = cauchywrappedpdf_width(f0,s);
	catch
	w=NaN;
	end
end % switch distribution
	% the density of thresholded patterns has no analytic expression,
	% so we estimate them as the average density over m patterns
	% generate random patterns
	b = real(ifft(sqrt(S).*fft(e)));
	% threshold
	b = b>mean(b,'all');
	s2 = var(b);
	% periodogram TODO normalize
	Shat = Lx./(s2*nx^2).*abs(fft(b-mean(b,'all'))).^2;
	% density estimate
	Sbar = mean(Shat,2);

	% estimate regularity based on various properties
	[reg(idx,:),leg_C,out(idx,jdx)] = regularity_measure(fx,S,distribution);
	% repeat estimating for density of thresholded patterns
	[reg_thresh(idx,:),leg_C,out_thresh(idx,jdx)] = regularity_measure(fx,Sbar,distribution);
end % for idx

splitfigure([2,3],[1,jdx],fflag);
cla();
loglog(Sxpc,reg,'linewidth',1);
if (jdx==1)
legend(leg_C,'Location','SouthEast')
end % if jdx==1
ylabel('Estimated regularity','interpreter','latex');
xlabel('Regularity $S_{xc}^+/\lambda_c$','interpreter','latex')
if (~pflag)
	title([upper(distribution(1)),distribution(2:end)])
end
hold on;
plot([1e-1,1e1],[1e-1,1e1],'k--') 

% thresholded estimates
splitfigure([2,3],[1,3+jdx],fflag);
cla();
loglog(Sxpc,reg_thresh,'linewidth',1);
if (jdx==1)
%legend(leg_C,'Location','SouthEast')
end % if jdx==1
ylabel('Estimated regularity','interpreter','latex');
xlabel('Regularity $S_{xc}^+/\lambda_c$','interpreter','latex')
if (~pflag)
	title([upper(distribution(1)),distribution(2:end)])
end
hold on;
plot([1e-1,1e1],[1e-1,1e1],'k--') 


end % for jdx
if (pflag)
	ps = 3;
	aspect=1.2;
	pdfprint(11,'img/regularity-measure-limit-gauss.pdf',ps,aspect);
	pdfprint(12,'img/regularity-measure-limit-laplace.pdf',ps,aspect);
	pdfprint(13,'img/regularity-measure-limit-cauchy.pdf',ps,aspect);
	pdfprint(14,'img/regularity-measure-limit-gauss-thresh.pdf',ps,aspect);
	pdfprint(15,'img/regularity-measure-limit-laplace-thresh.pdf',ps,aspect);
	pdfprint(16,'img/regularity-measure-limit-cauchy-thresh.pdf',ps,aspect);
end % if pflag

