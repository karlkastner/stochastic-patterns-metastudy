% Mon  1 Jul 10:39:45 CEST 2024
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
	pflag =0;
end
fflag = pflag;
sph  = 2; % 1.7

% characteristic frequency
fc = 1;
% density maximum
Sxpc = 1;

% spatial extent
L  = logspace(log10(0.55),1,50);
L  = logspace(0,log10(20),50);

% spatial resolution
dx = fc/100;

distribution_C = {'gauss','laplace','cauchy'};

reg = [];
reg_ = [];
entropy =[];
sd= [];
Rc = [];
% for each spatial extent
for idx=1:length(L)
	% spatial axis
	x  = (0:dx:L(idx))';
	% spectral axis
	fx = fourier_axis(x);
	% spectral resolution
	dfx = fx(2)-fx(1);
	% for each distribution
	for jdx=1:length(distribution_C)
	distribution = distribution_C{jdx};
	switch (distribution)
	case {'gauss'}
		% density parameter
		[f0,s]=normalwrappedpdf_mode2par(fc,Sxpc);
		% density
		Sx     = 0.5*normalwrappedpdf(fx,f0,s);
	case {'laplace'}
		[f0,s]=laplacewrappedpdf_mode2par(fc,Sxpc);
		Sx  = 0.5*laplacewrappedpdf(fx,f0,s);
	case {'cauchy'}
		[f0,s]=cauchywrappedpdf_mode2par(fc,Sxpc);
		Sx  = 0.5*cauchywrappedpdf(fx,f0,s);
	end

	% estimate regularity from various parameters
	[reg(idx,:,jdx), leg_C, out] = regularity_measure(fx,Sx,distribution);

	Rc(idx,jdx) = out.Rc;
	entropy(idx,jdx) = out.entropy;
	sd(idx,jdx) = out.sd;
if (0)
	if (1==jdx)
		figure(1)
		subplot(5,10,idx)
		cla
		plot(fx,S)
		hold on
		plot(fc_,Sxpc_,'*')
	end
end
	end
end
xlabel_str = 'Spatial extent L/\lambda_c'

xl_ = limits(dx);
figure(2);
clf
%yyaxis left
subplot(sph,4,1);
%plot(L,reg,'linewidth',1);
ylim([0,2.5]);
xlabel(xlabel_str);
%Spatial extent L/\lambda_c')
%yyaxis right
%ylim([0,1.8]);
ylabel('Estimated regularity S_c/\lambda_c');
%xlim(xl_);
legend('Gauss','Laplace','Cauchy');
%hline(reg(end,:))

subplot(sph,4,2);
plot(L,sd,'linewidth',1);
ylim([0,2.5]);
xlabel(xlabel_str);
ylabel('Standard deviation')
%hline(sd(end,:))
%xlim(xl_);

subplot(sph,4,3);
plot(L,entropy,'linewidth',1);
ylim([0,2.5]);
ylabel('Entropy');
xlabel(xlabel_str);
%hline(h(end,:))
%xlim(xl_);

subplot(sph,4,4);
plot(L,Rc);

for jdx=1:3;
 splitfigure([2,4],[1,4+jdx],fflag);
 plot(L,reg(:,:,jdx));
 xlabel(xlabel_str);
 title([upper(distribution_C{jdx}(1)), distribution_C{jdx}(2:end)]);
 ylim([0.5,1.5]);
 ylabel('Regularity estimate');
if (jdx == 1)
 legend(leg_C,'location','southeast');
end
end

if (pflag)
	ps = 3;
	aspect=1.2;
	pdfprint(15,'img/regularity-measure-spatial-extent-gauss.pdf',ps,aspect);
	pdfprint(16,'img/regularity-measure-spatial-extent-laplace.pdf',ps,aspect);
	pdfprint(17,'img/regularity-measure-spatial-extent-cauchy.pdf',ps,aspect);
end % if pflag

