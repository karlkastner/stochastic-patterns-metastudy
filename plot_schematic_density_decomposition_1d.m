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
%% plot synthetic periodic and stochastic patterns
%
function s = plot_schematic_density_decomposition_1d(meta)

	if (nargin()<1)
		meta = pattern_metastudy_metadata();
	end
	pflag = meta.pflag;
	fflag = pflag;
	issquare = true;
	
	cmap = colormap_kop();
	
	ppattern = 0.66;
	distribution = 'normal';
	
	% spatial extent
	L  = 20*[1,1];
	Lp = 10;
	% characteristic frequency, fc = 1/lambda_c
	fc = 1;
	% spatial resolution
	dx = 1./(fc*40);
	% number of points
	n  = round(L/dx);

	% maximum of the spectral density
	Sxpc = [0.5,1,2];
	Syc  = Sxpc;
	sy = 1.25;
	ns = 100;

	% coordinate axes in real space
	x = innerspace(0,L(1),n(1))-L(1)/2;
	y = innerspace(0,L(2),n(2))-L(2)/2;
	
	df = 1./L;

	% coordinate axes in frequency space
	[fx,fy,frr,tt] = fourier_axis_2d(L,n);
	%[fx,fy,frr]    = fourier_axis_2d(L,n);
		
	s=struct();
	% for isotropic, anisotropic
	for isiso=0:1
	% for periodic, periodic with noise, stochastic
	for idx=1:3
		% reset random number generator for exact reproducibility of figures
		rng(0);

	switch (idx)
	case {10}
		if (isiso)
			% periodic hexagonal pattern
			p = 1;
			q = 1;
			scale = false;
			a0 = 0;
			sbm = [];
			[b,x,y,Lx,Ly] = generate_isotropic_pattern(fc,n(1),L(1),a0,scale,sbm,p,q);
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
	case {20}
		% add noise to periodic patterns
		e=randn(n);
		b = ppattern*b./rms(b(:))+(1-ppattern)*e/rms(e(:));
		S = S/mean(S(:)) + ones(n);
	case {1,2,3}
		% stochastic pattern
		e  = randn(n);
		if (~isiso)
			% striped pattern
			switch (distribution)
			case {'normal'}
				[a,b] = normalmirroredpdf_mode2par(fc,0.5*Sxpc(idx));
				Sx    = normalmirroredpdf(fx,a,b);

				[a,b] = normpdf_mode2par(0,Syc(idx));
				Sy    = normpdf(fy,a,b);
			case {'lognormal'}
				[a,b] = lognpdf_mode2par(fc,Sxpc(idx));
				Sx    = lognmirroredpdf(fx,a,b);
				[a,b] = normpdf_mode2par(0,Syc(idx));
				Sy    = normpdf(fy,a,b);
			case {'gamma'}
				[a,b] = gampdf_mode2par(fc,Sxpc(idx));
				Sx    = gammirroredpdf(fx,a,b);
				%[a,b] = gamma_mode2par(1e-3,Sc(idx));
				[a,b]  = laplacepdf_max2par(Syc(idx))
				%Sy    = gampdf(abs(fy),a,b);
				Sy    = laplacepdf(fy,a,b);
			end
			Sxp = Sx.*(fx>=0);
			Sxy = cvec(Sx)*rvec(Sy);
		else
			% isotropic pattern
			switch (distribution)
			case {'normal'}
				[a,b] = normalmirroredpdf_mode2par(fc,0.5*Sxpc(idx));
				Sr_xy    = normalfoldedpdf(frr,a,b);
				Stc(idx) = Syc(idx)*2/pi;
				ct    = misesn_max2par(0.5*Stc(idx),6);
				St_xy = misesnpdf(tt,0,ct,6);
			case {'lognormal'}
				[a,b] = lognpdf_mode2par(fc,Sxpc(idx));
				Sr_xy    = lognpdf(frr,a,b);
				Stc(idx) = Syc(idx)*2/pi;
				ct  = misesn_max2par(0.5*Stc(idx),6);
				St_xy = misesnpdf(tt,0,ct,6);
			case {'gamma'}
				[a,b]    = gampdf_mode2par(fc,Sxpc(idx));
				Sr_xy    = gampdf(frr,a,b);
				Stc(idx) = Syc(idx)*2/pi;
				ct       = misesn_max2par(0.5*Stc(idx),6);
				St_xy    = misesnpdf(tt,0,ct,6);
			end
			Sxy = Sr_xy.*St_xy;
		end
		% transfer function
		T = sqrt(Sxy);
		% white (uncorrelated noise)
		bwhite = e;
		fwhite = ifft2(bwhite);
		% pattern
		b = real(ifft2(T.*fwhite));
	%	b = b/rms(b(:));
	case {4}
	end
	hatSxy = abs(fft2(b-mean(b(:)))).^2;
	% normalize
	hatSxy = hatSxy./(sum(hatSxy,'all')*df(1)*df(2));

	f_50 = fc;
	dfr = hypot(df(1),df(2));
	nf_test = round(0.25*f_50/dfr);
	fmsk = true(size(b));
	bmsk = true(size(b));
	fmsk = (frr<4*fc);
	bmsk = [];
	% TODO use spatial pattern analyis here
        [isperiodic, p_periodic, stati, out] = periodogram_test_periodicity_2d(...
					b, L, nf_test, bmsk, fmsk, ns);

	% autocorrelation
	Rxy = real(ifft2(Sxy));
	Rxy = Rxy/Rxy(1,1);
	Rx = mean(Rxy,2);
	Rx = Rx/Rx(1);
	Ry = mean(Rxy,1);
	Ry = Ry/Ry(1);
	[Rr,r]  = autocorr_radial(Rxy,L);
	[Rt,xt] = autocorr_angular(Rxy,L,101);

if (isiso)
	% the compoents have already been computed for anisotropic patterns

	% density components
	Sx = sum(Sxy,2)*df(2);
	% normalize
	Sx = Sx/(sum(Sx)*df(1));
	Sy  = sum(Sxy,1)*df(1);
	% normalize
	Sy  = Sy/(sum(Sy)*df(2));
end

	[Sr,fr] = periodogram_radial(Sxy,L);
	Sr      = Sr.normalized;
	[St,ft] = periodogram_angular(Sxy,L);
 
	% density along positive half-axis
	Sxp = 2*Sx.*(fx>0);
	Stp = 2*St.*(ft>-pi/2 & ft < +pi/2);

	if (isiso)
		[Sc_,mdx] = max(Sr);
		lc_ = 1./fr(mdx);
		reg = Sc_./lc_;
	else
		[Sxpc_,mdx] = max(Sxp);
		lc_ = 1./fx(mdx);
		reg = Sxpc_./lc_;
		
	end
	fprintf('iso %d p=%0.2f Sc/lc %0.2f\n',isiso,p_periodic,reg);
%		max(ratio(:)) 

	s(isiso+1,idx).x  = x;
	s(isiso+1,idx).y  = y;
	s(isiso+1,idx).b  = b;
	s(isiso+1,idx).Sxy = Sxy;
	s(isiso+1,idx).hatSxy = hatSxy;
	s(isiso+1,idx).r  = r;
	s(isiso+1,idx).Rr = Rr;
	s(isiso+1,idx).fx = fx;
	s(isiso+1,idx).Sx = Sx;
	s(isiso+1,idx).Sxp = Sxp;
	s(isiso+1,idx).fy = fy;
	s(isiso+1,idx).Sy = Sy;
	s(isiso+1,idx).fr = fr;
	s(isiso+1,idx).ft = ft;
	s(isiso+1,idx).Sr = Sr;
	s(isiso+1,idx).St = St;
	s(isiso+1,idx).Stp = Stp;

	end %  for idx

	% density on primary axis
	splitfigure([2,3],[3,1+3*isiso],fflag);
	cla;
	ls_C = {'-','-','-'};
	lw  = [1,1,1];
	for idx=1:3;
		if (isiso)
			plot(s(isiso+1,idx).fr,s(isiso+1,idx).Sr,ls_C{idx},'linewidth',lw(idx));
			lh =legend(num2str(cvec(Sxpc)));
			title(lh,'$S_{rc}/\lambda_c$','interpreter','latex');
		else
			plot(fftshift(s(isiso+1,idx).fx),fftshift(s(isiso+1,idx).Sxp),ls_C{idx},'linewidth',lw(idx));
			lh =legend(num2str(cvec(Sxpc)));
			title(lh,'$S_{xc}^+/\lambda_c$','interpreter','latex');
		end
		hold on;
	end
	 xlim([0,2.5]);
	 set(gca,'colororder',cmap)
	if (isiso)
	 xlabel('Wavenumber $k_r/k_c$','interpreter','latex');
	 ylabel('Density $S_r/k_c$','interpreter','latex');
	else
	 xlabel('Wavenumber $k_x/k_c$','interpreter','latex');
	 ylabel('Density $S_x^+/k_c$','interpreter','latex');
	if (issquare)
		axis square
	end	


	end

	% density along secondary axis	
	splitfigure([2,3],[3,2+3*isiso],fflag);
	cla();
	ls_C = {'-','-','-'};
	lw  = [1,1,1];
	for idx=1:3;
		if (isiso)
			plot(s(isiso+1,idx).ft,s(isiso+1,idx).Stp,ls_C{idx},'linewidth',lw(idx));
			lh =legend(arrayfun(@(x) sprintf('%g/\\pi',x),cvec(Stc*pi),'uniformoutput',false));
			title(lh,'$S_{{\theta}c}^+$','interpreter','latex');
		else
			plot(fftshift(s(isiso+1,idx).fy),fftshift(s(isiso+1,idx).Sy),ls_C{idx},'linewidth',lw(idx));
			lh =legend(num2str(cvec(Syc)));
			title(lh,'$S_{yc}/\lambda_c$','interpreter','latex');
		iSy = sum(s(isiso+1,idx).Sy)*(s(isiso+1,idx).fy(2)-s(isiso+1,idx).fy(1))

		end
		hold on;
	end
	 set(gca,'colororder',cmap)
	if (isiso)
	 xlabel('Angle $\theta$','interpreter','latex');
	 ylabel('Density $S_\theta^+$','interpreter','latex');
	 xlim([-pi,pi]/2)
	else
	 xlim([-2,2]);
	 xlabel('Wavenumber $k_y/k_c$','interpreter','latex');
	 ylabel('Density $S_y/k_c$','interpreter','latex');
	if (issquare)
		axis square
	end	
	end
	if (0)
	
	splitfigure([2,3],[3,3+3*isiso],fflag);
	cla();
	for idx=1:3;
		plot(NaN,NaN,ls_C{idx},'linewidth',lw(idx));
		hold on;
	end
	axis off
	set(gca,'colororder',cmap)
	legend('Periodic','Periodic + Noise','Stochastic')
	
	 splitfigure([2,3],[3,2+3*isiso],fflag);
	 cla();
	 for idx=1:3;
	 plot(s(isiso+1,idx).x1,s(isiso+1,idx).R1,ls_C{idx},'linewidth',lw(idx));
	 hold on;
	 end % for idx
	 xlim([0,2.5])
	 set(gca,'colororder',cmap)
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
		if (issquare)
			ps_ = 4;
		else
			ps = 3.5;
		end
		
%		pdfprint(2001,'img/pattern-decomposition-Sx.pdf',ps);
%		pdfprint(2002,'img/pattern-decomposition-Sy.pdf',ps);
%		pdfprint(2003,'img/pattern-decomposition-Sr.pdf',ps);
%		pdfprint(2004,'img/pattern-decomposition-St.pdf',ps);
%	
%		
%		pdfprint(101,'img/pattern-aniso-2d-periodic.pdf',ps);
%		pdfprint(107+2,'img/pattern-aniso-2d-periodic-with-noise.pdf',ps);
%		pdfprint(113+4,'img/pattern-aniso-2d-stochastic.pdf',ps);
%		
%		pdfprint(201,'img/pattern-iso-2d-periodic.pdf',ps);
%		pdfprint(207+2,'img/pattern-iso-2d-periodic-with-noise.pdf',ps);
%		pdfprint(213+4,'img/pattern-iso-2d-stochastic.pdf',ps);
%		
%		pdfprint(104,'img/periodogram-aniso-2d-periodic.pdf',ps);
%		pdfprint(110+2,'img/periodogram-aniso-2d-periodic-with-noise.pdf',ps);
%		pdfprint(116+4,'img/periodogram-aniso-2d-stochastic.pdf',ps);
%		pdfprint(105,'img/autocorrelation-aniso-2d-periodic.pdf',ps);
%		pdfprint(111+2,'img/autocorrelation-aniso-2d-periodic-with-noise.pdf',ps);
%		pdfprint(117+4,'img/autocorrelation-aniso-2d-stochastic.pdf',ps);
%		
%		pdfprint(204,'img/periodogram-iso-2d-periodic.pdf',ps);
%		pdfprint(210+2,'img/periodogram-iso-2d-periodic-with-noise.pdf',ps);
%		pdfprint(216+4,'img/periodogram-iso-2d-stochastic.pdf',ps);
%		pdfprint(205,'img/autocorrelation-iso-2d-periodic.pdf',ps);
%		pdfprint(211+2,'img/autocorrelation-iso-2d-periodic-with-noise.pdf',ps);
%		pdfprint(217+4,'img/autocorrelation-iso-2d-stochastic.pdf',ps);
%		
		pdfprint(31,'img/density-schematic-aniso-x.pdf',ps_);
		pdfprint(32,'img/schematic-density-aniso-y.pdf',ps_);
%		%pdfprint(33,'img/schematic-legend.pdf',ps);
		pdfprint(34,'img/density-schematic-iso-radial.pdf',ps);
		pdfprint(35,'img/schematic-density-iso-angular.pdf',ps);
%		
%		pdfprint(102,'img/pattern-1d-periodic.pdf',ps);
%		pdfprint(108,'img/pattern-1d-periodic-with-noise.pdf',ps);
%		pdfprint(114,'img/pattern-1d-stochastic.pdf',ps);
%		
%		pdfprint(103,'img/density-1d-periodic.pdf',ps);
%		pdfprint(109,'img/density-1d-periodic-with-noise.pdf',ps);
%		pdfprint(115,'img/density-1d-stochastic.pdf',ps);
%		
%		pdfprint(106,'img/autocorrelation-1d-periodic.pdf',ps);
%		pdfprint(112,'img/autocorrelation-1d-periodic-with-noise.pdf',ps);
%		pdfprint(118,'img/autocorrelation-1d-stochastic.pdf',ps);
	end
end

