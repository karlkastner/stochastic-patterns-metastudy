% Thu  4 Jul 20:57:06 CEST 2024
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
% estimate error and its components, i.e. bias and standard deviation
% for regularitie estimates based on the height of the mode
% for patterns with normal and log-normal spectral density
% spatial extent is kept constant

if (~exist('pflag','var'))
	pflag = 0;
end
fflag = pflag;

% characteristic wavelength
lc = 1;
% regularity
reg_x_r = logspace(-1,1,20)'
% density maximum
Sxpc_r = reg_x_r.*lc;
% characteristic frequency
fc = 1/lc;
% spatial extent
Lx = 20*fc;
Ly = Lx;
% spatial resolution
dx = fc/20;
dy = dx;
% spectral resolution
dfx = 1./Lx;
dfy = 1./Ly;
% number of grid points
nx = Lx/dx;
ny = Ly/dy;
% index of density maximum
idc = 1+Lx/lc;

% spatial axis
x = (0:nx-1)'*Lx/nx;
y = (0:ny-1)'*Ly/ny;
% spectral axis
fx = fourier_axis(x);
fy = fourier_axis(y);

% number of monte-carlo repetitions
m    = 400;

if (~exist('stat','var'))
	stat = struct();
	type_C = {'normal','lognormal'}
	
	for tdx=1:2
		%:length(type_C)
		type = type_C{tdx};
		% reset random number generator for reproducibility
		rng(0);
		
		se_reg_r = [];
		
		for rdx=1:length(reg_x_r)
		
		% regularity
		reg_x = reg_x_r(rdx);
		reg_y = reg_x;
	
		% density maximum	
		Sxpc = reg_x*lc;
		Sxy = reg_y*lc;
	
		%  distribution		
		switch (type)
		case {'normal'}
			% parameter
			[f0, sx] = normalwrappedpdf_mode2par(fc,Sxpc);
			% density
			Sx = 0.5*normalwrappedpdf(fx,f0,sx);
			Sxp = Sx.*(fx>=0);
			
			% cumulative density
			fdx = fx>=0;
			fx_=fx(fx>=0);
			C = cumsum(Sxp);
			C = C/C(end);
			fdx_=[true;C(2:end)~=C(1:end-1)];
			% quantiles
			f_me0 = interp1(C(fdx_),fx_(fdx_)+0.5./Lx,0.5,'linear')
			qr0 = diff(interp1(C(fdx_),fx_(fdx_)+0.5./Lx,[0.75,0.25],'linear'));
		
		case {'lognormal'}
			[a, b] = lognpdf_mode2par(fc,Sxpc);
			sx = lognpdf_std(a,b);
			Sx = 0.5*lognwrappedpdf(fx,a,b);
			f_me0 = logninv(0.5,a,b);
		end % switch

		% density parameter
		[~, sy] = normpdf_mode2par(0,Sxy);
		% density on y-axis
		Sy = normpdf(fy,0,sy);
		% two-dimensional density
		Sxy = cvec(Sx)*rvec(Sy);
		
		
		Sxyc = [];
		barSxpc = [];
		barSyc = [];
		sum_Sy2_df_ = [];
		f_me = [];
		% monte carlo repetision
		for idx=1:m
			% generate a random pattern
			e    = randn(nx,ny);
			Shatxy = Sxy.*abs(fft2(e)).^2/(nx*ny);
			% estimate the density
			barSx = sum(Shatxy,2)*dfy;
			% density on positive haf-axis
			barSxp = barSx.*(fx>=0);
			% normalize
			barSxp = barSxp/(sum(barSxp)*dfx);

			% TODO, use extreme3
			[barSxpc(idx,1),idc_] = max(barSxp);
			fc_(idx,1) = abs(fx(idc_));
			barSxpc(idx,2) = barSxp(idc);
			sum_Sy2_df_(idx,1) = sum(Shatxy(idc,:),2)*dfy;
			fdx = fx>=0;
			fx_=fx(fx>=0);
			C = cumsum(barSx(fdx));
			C = C/C(end);
			fdx_=[true;C(2:end)~=C(1:end-1)];
			f_me(idx,1) = interp1(C(fdx_),fx_(fdx_)+0.5./Lx,0.5,'linear');
			qr(idx,1) = diff(interp1(C(fdx_),fx_(fdx_)+0.5./Lx,[0.75,0.25],'linear'));
			%Sxyc(idx,1) = max(Shat,[],'all');
			%Sxyc(idx,2) = Shat(1,idc);	
		end % for idx
		reg_x_ = barSxpc.*fc_;
		% error of regularity estimate
		ereg   = (reg_x_ - Sxpc/fc);
		% error of characteristic frequency estimate
		efc    = fc_ - fc;
		% error of density maximum estimate
		eSxpc   = barSxpc - Sxpc; 
		
		% bias, standard deviation and standard error
		stat(tdx).b_qr_r(rdx,tdx)  = mean(qr - qr0);
		stat(tdx).sd_qr_r(rdx,tdx) = std(qr);
		stat(tdx).se_qr_r(rdx,tdx) = rms(qr - qr0);
		
		stat(tdx).b_f_me_r(rdx,tdx)    = mean(f_me - f_me0);
		stat(tdx).sd_f_me_r(rdx,tdx)    = std(f_me);
		stat(tdx).se_f_me_r(rdx,tdx)    = rms(f_me - f_me0);
		
		stat(tdx).sd_f_me_r(rdx,tdx) = std(f_me);
		stat(tdx).maxSxy_rat(rdx,tdx) = max(Sxy,[],'all');
		stat(tdx).sumSxy(rdx,tdx) = sum(Sxy,'all')*dfy*dfx;
		stat(tdx).bias_reg_r(rdx,1:2) = mean(ereg);
		stat(tdx).bias_Sxpc_r(rdx,1:2) = mean(eSxpc);
		stat(tdx).bias_fc_r(rdx,1)  = mean(efc);
		stat(tdx).sd_reg_r(rdx,1:2)   = std(eSxpc);
		stat(tdx).sd_Sxpc_r(rdx,1:2)   = std(eSxpc);
		stat(tdx).sd_fc_r(rdx,1)    = std(efc);
		stat(tdx).se_reg_r(rdx,1:2)   = rms(ereg);
		stat(tdx).se_Sxpc_r(rdx,1:2)   = rms(eSxpc);
		stat(tdx).se_fc_r(rdx,1)    = rms(efc);
		stat(tdx).corr_eSef(rdx,1)  = corr(efc,eSxpc(:,1));
		
		% test
		se_scx = sqrt(mean((barSxpc - Sxpc).^2));
		if (0)
		mu2 = mean(barSxpc.^2)
		mu = mean(barSxpc)
		se_scx = sqrt(mu2 - Sxpc^2)
		sd = std(barSxpc)
		%sd./mu
		sum_Sy2_df(1) = sum(Sy.^2)*dfy;
		sum_Sy2_df(2) = 1./sqrt(4*pi*sy^2)
		mu2_ = Sxpc^2*(sum_Sy2_df/Ly+1)
		se_scx = sqrt(mu2_ - Sxpc^2)
		se_scx = Sxpc*sqrt(sum_Sy2_df/Ly)
		se_scx = Sxpc*sqrt(1/(sqrt(4*pi*sy^2)*Ly))
		se_scx = Sxpc*sqrt(1/(sqrt(4*pi*1/(sqrt(2*pi)*reg_y*lc)^2)*Ly))
		se_scx = Sxpc*sqrt(1/(sqrt(2*1/(reg_y^2*lc^2)))*1/Ly)
		se_scx = Sxpc*sqrt(1/(sqrt(2)/(reg_y*lc))/Ly)
		end % if 0
		% analytic approximations
		stat(tdx).se_Sxpc_r(rdx,3) = reg_x*lc*sqrt(reg_y*lc/(sqrt(2)*Ly));
		stat(tdx).se_reg_r(rdx,3) = reg_x*sqrt(reg_y*lc/(sqrt(2)*Ly));
		%se_reg_r(rdx,:) = se_reg
		%se_fc_r(rdx,:) = rms(fc_-1./lc);
		%se = reg_x*reg_y*lc/Ly
		%se = Sxpc*sqrt(1/(sqrt(4*pi*(1/(reg_y*sqrt(2*pi)*lc).^2)*Ly)))
		%se = reg_x*sqrt(reg_x*(lc/Ly)^2)
		
		end % for rdx 
	end % for tdx 
end % if not exist stat
%end

for (tdx=1:length(stat))
type = type_C{tdx}

splitfigure([3,4],[tdx,1],fflag,'',100);
semilogx(reg_x_r,stat(tdx).se_reg_r./reg_x_r,'linewidth',1);
ylabel('Relative standard error $mse(reg_x)/reg_x$','interpreter','latex');
ylabel('$mse(reg_x)/reg_x$','interpreter','latex');

splitfigure([3,4],[tdx,2],fflag,'',100);
semilogx(reg_x_r,stat(tdx).se_Sxpc_r./Sxpc_r,'linewidth',1);
ylabel('Relative standard error $mse(S_{xc})/Sxpc$','interpreter','latex');
ylabel('$mse(S_{xc})/Sxpc$','interpreter','latex');

splitfigure([3,4],[tdx,3],fflag,'',100);
semilogx(reg_x_r,stat(tdx).se_fc_r/fc,'linewidth',1);
ylabel('Relative standard error $mse(f_c)/fc$','interpreter','latex');
ylabel('$mse(f_c)/fc$','interpreter','latex');

splitfigure([3,4],[tdx,5],fflag,'',100);
semilogx(reg_x_r,stat(tdx).sd_reg_r,'linewidth',1);
ylabel('Relative standard deviation $sd(reg_x)/reg_x$','interpreter','latex');
ylabel('$sd(reg_x)/reg_x$','interpreter','latex');

splitfigure([3,4],[tdx,6],fflag,'',100);
semilogx(reg_x_r,stat(tdx).sd_Sxpc_r,'linewidth',1);
ylabel('Relative standard deviation $sd(S_{xc})/S_{xc}$','interpreter','latex');
ylabel('$sd(S_{xc})/S_{xc}$','interpreter','latex');

splitfigure([3,4],[tdx,7],fflag,'',100);
semilogx(reg_x_r,stat(tdx).sd_fc_r,'linewidth',1);
ylabel('Relative standard deviation $sd(f_c)$','interpreter','latex');
ylabel('$sd(f_c)$','interpreter','latex');

splitfigure([3,4],[tdx,9],fflag,'',100);
semilogx(reg_x_r,stat(tdx).bias_reg_r./reg_x_r,'linewidth',1);
ylabel('Relative bias $bias(reg_x)/reg_x$','interpreter','latex');
ylabel('$bias(reg_x)/reg_x$','interpreter','latex');

splitfigure([3,4],[tdx,10],fflag,'',100);
semilogx(reg_x_r,stat(tdx).bias_Sxpc_r./Sxpc_r,'linewidth',1);
ylabel('Relative Bias $bias(S_{xc})/S_{xc}$','interpreter','latex');
ylabel('$bias(S_{xc})/S_{xc}$','interpreter','latex');

splitfigure([3,4],[tdx,11],fflag,'',100);
semilogx(reg_x_r,stat(tdx).bias_fc_r./fc,'linewidth',1);
ylabel('Relative Bias $bias(f_c)/f_c$','interpreter','latex');
ylabel('$bias(f_c)/f_c$','interpreter','latex');

%semilogx(reg_x_r,stat(tdx).se_reg_r./reg_x_r,'linewidth',1);
%xlabel('Regularity Sc/lc');
%ylabel('Relative standard error se_{reg_x}/reg_x');
%legend('approximation','MC-estimate')

splitfigure([3,4],[tdx,4],fflag,'',100);
semilogx(reg_x_r,stat(tdx).se_reg_r./reg_x_r./sqrt(1/sqrt(2)*reg_x_r*lc/Ly),'linewidth',1);
xlabel('Regularity S_{xc}/\lambda_c');
ylabel({'Normalized standard error', '$se_{reg_x} \sqrt{L_y/(\sqrt{2} \lambda_c reg_x^3)}$'},'interpreter','latex');
legend('MC, S_{xc} and f_c estimated','MC, S_{xc} estimated, f_c exact','limit reg_x \rightarrow \infty')

splitfigure([3,4],[tdx,8],fflag,'',100);
semilogx(reg_x_r,stat(tdx).sd_reg_r./reg_x_r./sqrt(1/sqrt(2)*reg_x_r*lc/Ly),'linewidth',1);
xlabel('Regularity S_{xc}/\lambda_c');
ylabel({'Normalized standard error', '$se_{reg_x} \sqrt{L_y/(\sqrt{2} \lambda_c reg_x^3)}$'},'interpreter','latex');
%legend('MC-estimate','MC-estimate, fc exact','limit reg_x \rightarrow \inf')

splitfigure([3,4],[tdx,12],fflag,'',100);
semilogx(reg_x_r,stat(tdx).bias_reg_r./reg_x_r./sqrt(1/sqrt(2)*reg_x_r*lc/Ly),'linewidth',1);
xlabel('Regularity S_{xc}/\lambda_c');
ylabel({'Normalized standard error', '$se_{reg_x} \sqrt{L_y/(\sqrt{2} \lambda_c reg_x^3)}$'},'interpreter','latex');
%legend('MC-estimate','MC-estimate, fc exact','limit reg_x \rightarrow \inf')


splitfigure([3,4],[10+tdx,1],fflag,'',100);
semilogx(reg_x_r,stat(tdx).corr_eSef,'linewidth',1);
ylabel('Correlation $ef_c, eSxpc$','interpreter','latex');

splitfigure([3,4],[10+tdx,2],fflag,'',100);
semilogx(reg_x_r,stat(tdx).se_f_me_r,'linewidth',1);

splitfigure([3,4],[10+tdx,2+4],fflag,'',100);
semilogx(reg_x_r,stat(tdx).sd_f_me_r,'linewidth',1);
%ylabel('Correlation $ef_c, eSxpc$','interpreter','latex')

splitfigure([3,4],[10+tdx,2+8],fflag,'',100);
semilogx(reg_x_r,stat(tdx).b_f_me_r,'linewidth',1);

%ylabel('Correlation $ef_c, eSxpc$','interpreter','latex')
splitfigure([3,4],[10+tdx,3],fflag,'',100);
semilogx(reg_x_r,stat(tdx).se_qr_r,'linewidth',1);
ylabel('se_qr')

splitfigure([3,4],[10+tdx,3+4],fflag,'',100);
semilogx(reg_x_r,stat(tdx).sd_qr_r,'linewidth',1);
ylabel('sd_qr')

splitfigure([3,4],[10+tdx,3+8],fflag,'',100);
semilogx(reg_x_r,stat(tdx).b_qr_r,'linewidth',1);
ylabel('b_qr')

%subplot(2,2,3);
%$splitfigure([2,3],[1,4],fflag,'',100)
%semilogx(reg_x_r,bias_r./reg_x_r,'linewidth',1)

%splitfigure([2,3],[1,5],fflag,'',100)
%semilogx(reg_x_r,se_fc_r./fc,'linewidth',1)

tdx
if(1)
if (~pflag)
	pdfprint(100*tdx,['img/uncertainty-regularity-estimate-overview-',type,'.pdf'],1);
else
	ps = 3.5;
	pdfprint(100*tdx+1,['img/uncertainty-estimate-',type,'-mse-regx.pdf'],ps);
	pdfprint(100*tdx+2,['img/uncertainty-estimate-',type,'-mse-Sxpc.pdf'],ps);
	pdfprint(100*tdx+3,['img/uncertainty-estimate-',type,'-mse-fc.pdf'],ps);
	pdfprint(100*tdx+4,['img/uncertainty-estimate-',type,'-convergence.pdf'],ps);
	pdfprint(100*tdx+1+4,['img/uncertainty-estimate-',type,'-sd-regx.pdf'],ps);
	pdfprint(100*tdx+2+4,['img/uncertainty-estimate-',type,'-sd-Sxpc.pdf'],ps);
	pdfprint(100*tdx+3+4,['img/uncertainty-estimate-',type,'-sd-fc.pdf'],ps);
	pdfprint(100*tdx+1+8,['img/uncertainty-estimate-',type,'-bias-regx.pdf'],ps);
	pdfprint(100*tdx+2+8,['img/uncertainty-estimate-',type,'-bias-Sxpc.pdf'],ps);
	pdfprint(100*tdx+3+8,['img/uncertainty-estimate-',type,'-bias-fc.pdf'],ps);
end
end

end % for tdx 


