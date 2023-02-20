% 2022-12-09 10:51:33.398412971 +0100
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
% plot synthetic periodic and stochastic patterns
%
function pattern_schematic_plot(meta)

	if (nargin()<1)
		meta = pattern_metadata();
	end
	pflag = meta.pflag;
	fflag = pflag;
	
	ppattern = 0.66;
	islogn   = false;
	
	% spatial extent
	L  = 20*[1,1];
	Lp = 10;
	% characteristic frequency, fc = 1/lambda_c
	fc = 1;
	% spatial resolution
	dx = 1./(fc*20);
	% number of points
	n  = round(L/dx);

	% maximum of the spectral density
	Sc = 1.5;
	sy = 1.25;
	ns = 100;
		
	s=struct();
	% for isotropic, anisotropic
	for isiso=0:1
		% coordinate axes in real space
		x = innerspace(0,L(1),n(1))-L(1)/2;
		y = innerspace(0,L(2),n(2))-L(2)/2;
		% coordinate axes in frequency space
		[fx,fy,frr] = fourier_axis_2d(L,n);
	% for periodic, periodic with noise, stochastic
	for idx=1:3
		% reset random number generator for exact reproducibility of figures
		rng(0);

	switch (idx)
	case {1}
		if (isiso)
			% periodic hexagonal pattern
			p = 1;
			q = 1;
			scale = false;
			a0 = 0;
			sbm = [];
			[b,x,y,Lx,Ly] = hexagonal_pattern(fc,n(1),L(1),a0,scale,sbm,p,q);
			x=x-mean(x);
			y=y-mean(y);
			L = [Lx,Ly];
			n = size(b);
			[fx,fy,frr] = fourier_axis_2d(L,n);
		else
			% periodic striped pattern
		b = sin(2*pi*x*fc + 0.*y')';
		end
		S = abs(fft2(b)).^2;
		S(1,1) = 0;
	case {2}
		% add noise to periodic patterns
		e=randn(n);
		b = ppattern*b./rms(b(:))+(1-ppattern)*e/rms(e(:));
		S = S/mean(S(:)) + ones(n);
	case {3}
		% stochastic pattern
		e  = randn(n);
		if (~isiso)
			% striped pattern
			if (islogn)
				[a,b] = logn_mode2param(fc,Sc);
				Sx    = lognpdf(fx,a,b);
				Sy    = normpdf(fy,0,0.5./fc);
			else
				[a,b] = gamma_mode2par(fc,Sc);
				Sx    = gampdf(fx,a,b);
				Sy    = gampdf(abs(fy),1,0.25);
			end
			S = cvec(Sx)*rvec(Sy);
		else
			% isotropic pattern
			if (islogn)
				[a,b] = logn_mode2param(fc,Sc);
				S = lognpdf(frr,a,b);
			else
				[a,b] = gamma_mode2par(fc,Sc);
				S = gampdf(frr,a,b);
			end
		end
		% transfer function
		T = sqrt(S);
		% white (uncorrelated noise)
		bwhite = e;
		fwhite = ifft2(bwhite);
		% pattern
		b = real(ifft2(T.*fwhite));
	%	b = b/rms(b(:));
	case {4}
	end
	[fx,fy,frr] = fourier_axis_2d(L,n);
	Shat = abs(fft2(b-mean(b(:)))).^2;
	df = 1./L;
	Shat = 2*Shat./(sum(Shat(:))*df(1)*df(2));

	f_50 = fc;
	dfr = df(1);
	nf_test = round(0.25*f_50/dfr);
	fmsk = true(size(b));
	bmsk = true(size(b));
	fmsk = (frr<4*fc);
	bmsk = [];
	% TODO use spatial pattern analyis here
        [p_periodic, stati, out] = periodogram_test_periodicity_2d(...
					b, nf_test, bmsk, fmsk, ns);
	
	R = real(ifft2(S));
	R = R/R(1,1);
	Sx = mean(S,2);
	fdx = (fx>=0);
	df = (fx(2)-fx(1));
	Sx = Sx/(sum(Sx(fdx))/L(1));

	Sy  = mean(S,1);
	fdx = true(size(fy));
	Sy  = Sy/(sum(Sy(fdx))/L(1));

	Rx = mean(R,2);
	Rx = Rx/Rx(1);
	Ry = mean(R,1);
	Ry = Ry/Ry(1);

	[Sr,fr] = periodogram_radial(S,L);
	Sr      = Sr.normalized;
	[Rr,r]  = autocorr_radial(R,L);
	[St,ft] = periodogram_angular(S,L);
	[Rt,xt] = autocorr_angular(R,L,101);

	if (isiso)
		[Sc_,mdx] = max(Sr);
		lc_ = 1./fr(mdx);
		reg = Sc_./lc_;
	else
		[Sc_,mdx] = max(Sx);
		lc_ = 1./fx(mdx);
		reg = Sc_./lc_;
		
	end
	fprintf('iso %d p=%0.2f Sc/lc %0.2f\n',isiso,p_periodic,reg);
%		max(ratio(:)) 

	s(isiso+1,idx).x  = x;
	s(isiso+1,idx).y  = y;
	s(isiso+1,idx).b  = b;
	s(isiso+1,idx).Shat = Shat;
	s(isiso+1,idx).r  = r;
	s(isiso+1,idx).Rr = Rr;
	s(isiso+1,idx).fx = fx;
	s(isiso+1,idx).Sx = Sx;
	s(isiso+1,idx).fy = fy;
	s(isiso+1,idx).Sy = Sy;
	s(isiso+1,idx).fr = fr;
	s(isiso+1,idx).ft = ft;
	s(isiso+1,idx).Sr = Sr;
	s(isiso+1,idx).St = St;
	
	% density on secondary axis
	splitfigure([3,8],[10+isiso*10,(idx-1)*8+7],fflag);
	cla
	if (isiso)
		plot(ft,St);
		ylabel('S_t');
		xlabel('\theta');
	else
		plot(fy,Sy);
		ylabel('S_y');
		xlabel('k_y/\k_c');
	end
	
	% autocorrelation
	splitfigure([3,8],[10+isiso*10,(idx-1)*8+8],fflag);
	cla
	if (isiso)
		plot(xt,Rt);
		ylabel('R_t');
		xlabel('\theta');
	else
		plot(y,Ry);
		ylabel('R_y');
		xlabel('y/\lambda_c');
	end

	% pattern (thresholded)	
	splitfigure([3,8],[10+isiso*10,(idx-1)*8+1],fflag);
	cla
	if (isiso)
		p_ = 0.8;
	else
		p_ = 0.55;
	end
	imagesc(x*fc,y*fc,b'>quantile(b(:),p_));
	axis([0,Lp,0,Lp]*(1+sqrt(eps)))
	axis square
	set(gca,'xtick',0:2:(Lp-1),'ytick',0:2:(Lp-1));
	xlabel('Distance $x/\lambda_c$','interpreter','latex')';
	ylabel('Distance $y/\lambda_c$','interpreter','latex')';
	colormap(flipud(gray))
	axis xy
	
	% pattern along x and y
	splitfigure([3,8],[10+isiso*10,(idx-1)*8+2],fflag);
	cla
	if (~isiso)
	plot(x*fc,b(:,1),'linewidth',1);
	else
	plot(y*fc,b(1,:),'linewidth',1);
	end
	ylim([-1.4,1.4*3])
	daspect([1,1,1]);
	set(gca,'ytick',[],'xtick',0:2:10,'xlim',[-sqrt(eps),10+sqrt(eps)]);
	xlabel('Distance $x/\lambda_c$','interpreter','latex')
	
	splitfigure([3,8],[10+isiso*10,(idx-1)*8+3],fflag);
	cla
	if (isiso)
		plot(fftshift(fr/fc),(fftshift(Sr))*fc,'linewidth',1);
	else
		plot(fftshift(fx/fc),(fftshift(Sx))*fc,'linewidth',1);
	end
	xlim([0,2.5]) 
	xlabel('Wavenumber $k_x/k_c$','interpreter','latex');
	ylabel('Density $S_x/k_c$','interpreter','latex');
	

	% periodogram
	splitfigure([3,8],[10+isiso*10,(idx-1)*8+4],fflag);
	cla
	
	if (idx==1 || idx == 2)
		Shat_ = Shat;
		flag = double(Shat_ < 0.9*max(Shat_(:)));
		Shat_ = (Shat_.*flag) + circfilt2(Shat_.*(1-flag),1.001);
		Shat_ = max(0,Shat_);
		Shat_ = log10(Shat_);
	else
		Shat_ = log10(Shat);
	end
	Shat_(1,1) = log10(eps);
	
	df = 1./L;
	scale = log10(L(1)*L(2));
	Shat = scale*Shat_;
	imagesc(fftshift(fx)/fc,fftshift(fy)/fc,fftshift(Shat_)');
	axis(2.5*[-1,1,-1,1]);
		
	if (idx ~= 3)
	if (isiso)
		caxis(max(Shat_(:))+[-4.05,0]);
	else
		caxis(max(Shat_(:))+[-4.9,0]);
	end
	else
		caxis(max(Shat_(:))+[-2,0]);
	end
	axis square
	
	colormap(flipud(gray))
	xlabel('Wavenumber $k_x/k_c$','interpreter','latex');
	ylabel('Wavenumber $k_x/k_c$','interpreter','latex');
	
	% 2d-autocorrelation
	splitfigure([3,8],[10+isiso*10,(idx-1)*8+5],fflag);
	cla
	imagesc(x,y,fftshift(R)')
	xlim([-2.5,2.5]);
	ylim([-2.5,2.5]);
	axis square
	xlabel('Lag Distance x / \lambda_c');
	ylabel('Lag Distance y / \lambda_c');
	
	% autocorrelation along primary axis
	splitfigure([3,8],[10+isiso*10,(idx-1)*8+6],fflag);
	cla
	if (isiso)
		plot(r*fc,Rr,'linewidth',1);
		ylim([-0.5,1]);
		ylabel('R_r');
		xlabel('Lag distance $r/\lambda_c$','interpreter','latex');
	else
		plot(x*fc,Rx,'linewidth',1);
		ylabel('R_x');
		ylim([-1,1]);
		xlabel('Lag distance $x/\lambda_c$','interpreter','latex');
	end % if isiso
	xlim([0,2.5])
	if (0)
		hold on
		R1_ = autocorr_radial_hexagonal_pattern(x1,1./fc);
		plot(x1*fc,R1_);
	end

	if (idx == 3)
		fflag = 1;
		splitfigure([2,2],[200,2*isiso+1],fflag);
		cla
		if (isiso)
		plot(fr/fc,Sr,'linewidth',1); hold on		
		plot(-fr/fc,Sr,'k--','linewidth',1);
		ylim([0,1.55]);
		xlabel('Wavenumber k_r/k_c');
		ylabel('Density S_r/\lambda_c');
%		plot
		else
		plot(fx/fc,Sx,'linewidth',1); hold on		
		plot(-fx/fc,Sx,'k--','linewidth',1);
		xlabel('Wavenumber k_x/k_c');
		ylabel('Density S_x/\lambda_c');
		ylim([0,1.55]);
		end
		xlim([-2.5,2.5]);
		splitfigure([2,2],[200,2*isiso+2],fflag);
		cla
		if (isiso)
		fdx = abs(ft)<= pi/2;
		plot(ft(fdx),St(fdx),'linewidth',1); hold on		
		fdx = (ft)>= pi/2;
		plot(ft(fdx),St(fdx),'k--','linewidth',1); hold on	
		fdx = (ft)<=- pi/2;
		plot(ft(fdx),St(fdx),'k--','linewidth',1); hold on	
		xlim([-pi,pi]);
		set(gca,'xtick',[-pi/2,0,pi/2],'xticklabel',{'-\pi/2','0','\pi/2'})
		ylim([0,1.55])	
		xlabel('Angle \theta');
		ylabel('Density S_\theta');
		daspect([2,1,1]); 
		camroll(90);
		else
		fdx = fy>=0;
		plot(fy(fdx)/fc,Sy(fdx),'linewidth',1); hold on		
		plot(-fy(fdx)/fc,Sy(fdx),'k--','linewidth',1);
		
		xlabel('Wavenumber k_y/k_c');
		ylabel('Density S_y/\lambda_c');
		xlim([-2.5,2.5])
		ylim([0,1.55])	
		camroll(90);
		end
	
		fflag=0;
	end % if idx == 3

	end %  for idx

	% density on primary axis
	splitfigure([2,3],[3,1+3*isiso],fflag);
	 cla;
	ls_C = {'-','--','-'};
	lw  = [1,1.5,1];
	for idx=1:3;
		if (isiso)
			plot(s(isiso+1,idx).fr,s(isiso+1,idx).Sr,ls_C{idx},'linewidth',lw(idx));
		else
			plot(s(isiso+1,idx).fx,s(isiso+1,idx).Sx,ls_C{idx},'linewidth',lw(idx));
		end
		hold on;
	end
	 xlim([0,2.5]);
	 set(gca,'colororder',meta.colormap)
	if (isiso)
	 xlabel('Wavenumber $k_r/k_c$','interpreter','latex');
	 ylabel('Density $S_r/k_c$','interpreter','latex');
	else
	 xlabel('Wavenumber $k_x/k_c$','interpreter','latex');
	 ylabel('Density $S_x/k_c$','interpreter','latex');
	end

	% density along secondary axis	
	splitfigure([2,3],[3,2+3*isiso],fflag);
	cla();
	ls_C = {'-','--','-'};
	lw  = [1,1.5,1];
	for idx=1:3;
		if (isiso)
			plot(s(isiso+1,idx).ft,s(isiso+1,idx).St,ls_C{idx},'linewidth',lw(idx));
		else
			plot(s(isiso+1,idx).fy,s(isiso+1,idx).Sy,ls_C{idx},'linewidth',lw(idx));
		end
		hold on;
	end
	 set(gca,'colororder',meta.colormap)
	if (isiso)
	 xlabel('Angle $\theta$','interpreter','latex');
	 ylabel('Density $S_\theta$','interpreter','latex');
	 xlim([-pi,pi]/2)
	else
	 xlim([0,2.5]);
	 xlabel('Wavenumber $k_y/k_c$','interpreter','latex');
	 ylabel('Density $S_y/k_c$','interpreter','latex');
	end
	if (0)
	
	splitfigure([2,3],[3,3+3*isiso],fflag);
	cla();
	for idx=1:3;
		plot(NaN,NaN,ls_C{idx},'linewidth',lw(idx));
		hold on;
	end
	axis off
	set(gca,'colororder',meta.colormap)
	legend('Periodic','Periodic + Noise','Stochastic')
	
	 splitfigure([2,3],[3,2+3*isiso],fflag);
	 cla();
	 for idx=1:3;
	 plot(s(isiso+1,idx).x1,s(isiso+1,idx).R1,ls_C{idx},'linewidth',lw(idx));
	 hold on;
	 end % for idx
	 xlim([0,2.5])
	 set(gca,'colororder',meta.colormap)
	 xlabel('Lag distance $x/\lambda_c$','interpreter','latex');
	 ylabel('Density $S_x/k_c$','interpreter','latex');
	if (isiso)
	 ylabel('Autocorrelation $R_r$','interpreter','latex');
	else
	 ylabel('Autocorrelation $R_x$','interpreter','latex');
	end
	
	end
	
	end
	
	if (pflag)
		ps = 3.5;
		
		pdfprint(2001,'img/pattern-decomposition-Sx.pdf',ps);
		pdfprint(2002,'img/pattern-decomposition-Sy.pdf',ps);
		pdfprint(2003,'img/pattern-decomposition-Sr.pdf',ps);
		pdfprint(2004,'img/pattern-decomposition-St.pdf',ps);
	
		
		pdfprint(101,'img/pattern-aniso-2d-periodic.pdf',ps);
		pdfprint(107+2,'img/pattern-aniso-2d-periodic-with-noise.pdf',ps);
		pdfprint(113+4,'img/pattern-aniso-2d-stochastic.pdf',ps);
		
		pdfprint(201,'img/pattern-iso-2d-periodic.pdf',ps);
		pdfprint(207+2,'img/pattern-iso-2d-periodic-with-noise.pdf',ps);
		pdfprint(213+4,'img/pattern-iso-2d-stochastic.pdf',ps);
		
		pdfprint(104,'img/periodogram-aniso-2d-periodic.pdf',ps);
		pdfprint(110+2,'img/periodogram-aniso-2d-periodic-with-noise.pdf',ps);
		pdfprint(116+4,'img/periodogram-aniso-2d-stochastic.pdf',ps);
		pdfprint(105,'img/autocorrelation-aniso-2d-periodic.pdf',ps);
		pdfprint(111+2,'img/autocorrelation-aniso-2d-periodic-with-noise.pdf',ps);
		pdfprint(117+4,'img/autocorrelation-aniso-2d-stochastic.pdf',ps);
		
		pdfprint(204,'img/periodogram-iso-2d-periodic.pdf',ps);
		pdfprint(210+2,'img/periodogram-iso-2d-periodic-with-noise.pdf',ps);
		pdfprint(216+4,'img/periodogram-iso-2d-stochastic.pdf',ps);
		pdfprint(205,'img/autocorrelation-iso-2d-periodic.pdf',ps);
		pdfprint(211+2,'img/autocorrelation-iso-2d-periodic-with-noise.pdf',ps);
		pdfprint(217+4,'img/autocorrelation-iso-2d-stochastic.pdf',ps);
		
		pdfprint(31,'img/density-schematic-aniso.pdf',ps);
		pdfprint(32,'img/schematic-density-aniso-y.pdf',ps);
		%pdfprint(33,'img/schematic-legend.pdf',ps);
		pdfprint(34,'img/density-schematic-iso.pdf',ps);
		pdfprint(35,'img/schematic-density-iso-angular.pdf',ps);
		
		if (0)
		pdfprint(102,'img/pattern-1d-periodic.pdf',ps);
		pdfprint(108,'img/pattern-1d-periodic-with-noise.pdf',ps);
		pdfprint(114,'img/pattern-1d-stochastic.pdf',ps);
		
		pdfprint(103,'img/density-1d-periodic.pdf',ps);
		pdfprint(109,'img/density-1d-periodic-with-noise.pdf',ps);
		pdfprint(115,'img/density-1d-stochastic.pdf',ps);
		
		pdfprint(106,'img/autocorrelation-1d-periodic.pdf',ps);
		pdfprint(112,'img/autocorrelation-1d-periodic-with-noise.pdf',ps);
		pdfprint(118,'img/autocorrelation-1d-stochastic.pdf',ps);
		end
	end
end

