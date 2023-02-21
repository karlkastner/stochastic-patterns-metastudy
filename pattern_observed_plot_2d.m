% 2021-12-03 10:51:02.490036669 +0100
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
% plot the remotely sensed patterns and their spectral densities displayed in figure 2
%
function sp_a = pattern_observed_plot_2d(meta)
	if (nargin()<1)
		meta = pattern_analysis_metadata();
	end
	ut = true;
	fflag = meta.pflag;
	pline = false;
	rotate = true;
	lineflag = 0;
	xscale = true;

	f_C = { 
		'img/this-publication/2a_11.33051_28.35545_model_0_mercator.png', 100, 2.5, ...
		'img/2d/2/pattern_2d_+11.53386_+027.92788.png',  60, 3.5,
	};


	sp_a = Spatial_Pattern();

	% for both the spotted and the banded pattern
	for jdx=1:size(f_C,1)
	rmax     = f_C{jdx,2};
	fxscale  = f_C{jdx,3};

	% set options
	sp = Spatial_Pattern();
	sp.opt.rmax = rmax;

	% read image
	sp.imread(f_C{jdx,1});
	if (jdx == 1)
		sp.L = 2*sp.n;
	end
	n   = size(sp.b);
	L   = sp.L;

	% for repeat the analysis twice set kdx=1:2
	% the first time for the orignal pattern,
	% the second time for an idealized periodic pattern with the same wavelength
	for kdx=1
	% generate idealized periodic patterns at second iterationg
	if (2 == kdx)
		if (1==jdx)
			[b] = band_pattern(fc(1),n,L,sp.stat.fct);
		else 
			[b] = hexagonal_pattern(fc(1),n,L,pi/6);
		end 
		sp.b = b;
	end
	
	% analyze the pattern
	rng(0)
	sp.analyze_grid();
	banded = ~sp.stat.isisotropic();
	%sp.fit_parametric_densities();
	f.r    = sp.f.r;
	f.x    = fftshift(sp.f.x);
	f.y    = fftshift(sp.f.y);
	fc(kdx) = sp.stat.fc.radial.clip;
	reg    = sp.stat.Sc.x.clip.*fc(kdx);

	fprintf('P-value of periodicity test: %f\n',sp.stat.p_periodic);
	fprintf('Wavelength: %f\n',1./fc(kdx));
	fprintf('Regularity: %f\n',reg(kdx));

	Shat = sp.S.hat;

	% plot the pattern
	[~,~,tp] = splitfigure([2,4],[jdx,1],fflag,'',10.^kdx,ut);
	cla();
	b_ = sp.b_square;
	if (xscale)
		if (sp.stat.isisotropic)
			lc = 1./sp.stat.fc.radial.clip;
		else
			lc = 1./sp.stat.fc.x.clip;
		end
		xs=1./lc;
	end
	sp.plot('b');

	if(~xscale)
	set(gca,'xtick',-1500:250:1500,'ytick',-1500:250:1500);
	end
	clim_ = quantile(b_(:),[0.2,0.8]);
	clim(clim_);
	hold on
	try
	if (lineflag)
	plot(y(round(yc_([1,end]))),x(round(xc_([1,end]))),'r','linewidth',1.5);
	end
	catch
	end
	L_ = min(sp.L);
	axis xy 
	axis tight
	xlim([0,1]*L_/lc);
	ylim([0,1]*L_/lc);
	axis square
	if (fflag)
		title('');
		colormap gray
		axis on;
		if (xscale)
			xlabel('$x/\lambda_c$','interpreter','latex');
			ylabel('$y/\lambda_c$','interpreter','latex');
		else
			xlabel('Easting (m)'); ylabel('Northing (m)');
		end
	else
	title('Pattern (biomass proxy)','interpreter','latex');
	end
if (0)
	ax2 = axes(tp);
	xs = 1e-3;
	imagesc(ax2,xs*(sp.x-sp.x(1)),xs*(sp.y-sp.y(1)),sp_b_','alphadata',1);
	ax2.XAxisLocation = 'top';
	ax2.YAxisLocation = 'right';
	ax2.Color = 'none';
	ax2.Box = 'off';
	axis(ax2,'tight')
	xlim(ax2,[0,L_]*xs);
	ylim(ax2,[0,L_]*xs);
	axis(ax2,'xy');
	xlabel('$x$/km','interpreter','latex');
	ylabel('$y$/km','interpreter','latex');
end
	% plot the 2D-periodogram
	splitfigure([2,4],[jdx,2],fflag,'',10.^kdx,ut,tp);
	cla();
	% rotate, to make bands parallel to y
	sp.plot('S.rot.clip');
	% note: surface not properly plotted as svg,'edgecolor','none');
	if (~meta.pflag) % && ~rotate)
		hold on
		%plot(-sp.stat.fc.y.clip/fc(kdx),-sp.stat.fc.x.clip/fc(kdx),'r*')
	end
	axis(fxscale*[-1,1,-1,1]);
	if (lineflag)
		hold on
		if (rotate)
			plot([0,0],[-3,3],'r--','linewidth',1.5);
			plot([0,0],[0,-3],'r','linewidth',1.5);
			plot([-3,3],-[1,1],'b--','linewidth',1.5);
			plot([0,3],-[1,1],'b','linewidth',1.5);
		else
			plot(f.y(round(jc([1,end])))/fc(kdx),f.x(round(ic([1,end])))/fc(kdx),'r','linewidth',1.5);
		end
	end
	if (fflag)
		axis on
	else
		axis off;
		title('Periodogram $\hat S$','interpreter','latex');
	end

	% plot the 2D-autocorrelation
	splitfigure([2,4],[jdx,3],fflag,'',10.^kdx,ut,tp);
	cla();
	sp.plot('R.clip');
	title('Autocorrelation $R$','interpreter','latex');
	axis(4*[-1,1,-1,1]);

	% plot the 2d spectral density
	splitfigure([2,4],[jdx,4],fflag,'',10.^kdx, ut,tp);
	cla();
	sp.plot('S.rot.bar');
	axis(fxscale*[-1,1,-1,1]);
	hold on 
	if (lineflag)
%		plot(f.y(round(jc_([1,end])))/fc(kdx),f.x(round(ic_([1,end])))/fc(kdx),'b','linewidth',1.5);
	end
	axis off
	if (fflag)
		%axis off
		title('');
	else
		title('Spectral Density $\hat S$','interpreter','latex');
		cl = colorbar();
		title(cl,'${\bar S}/\lambda_c^2$','interpreter','latex')
	end

	% plot pattern biomass proxy along the transect
if (0)
	splitfigure([2,4],[jdx,5],fflag,'',10.^kdx, ut,tp);
	cla();
	xl   = [0.5439, 3.8077];
	ylabel('Normalized biomas proxy')
	xlabel('$r/\lambda_c$','interpreter','latex');
	xlim([0,range(xl)]);
end

	% 1D autocorrelation in direction perpendicular to bands or radial for spotted pattern
	splitfigure([2,4],[jdx,6],fflag,'',10.^kdx, ut,tp);
	cla();
	if (sp.stat.isisotropic)
		plot(sp.x*fc(kdx),fftshift(sp.R.rot.x.clip));
		ylabel('$R_x$','rot',0,'interpreter','latex');
	else
		plot(sp.r*fc(kdx),sp.R.radial.clip);
		ylabel('$R_r$','rot',0,'interpreter','latex');
	end
	xlim([0,2.5]); %fc(kdx)*rmax]);
	xlabel('$x / \lambda_c$','interpreter','latex')
	
	% 1D density in direction perpendicular to bands or radial for spotted pattern
	splitfigure([2,4],[jdx,7],fflag,'',10.^kdx, ut,tp);
	cla();
	f.x = sp.f.x;
	fdx = f.x>=0;
	if (~sp.stat.isisotropic)
		sp.plot('S.rot.x.hat');
		hold on
		if (0)
		sp.plot('S.rot.x.brownian_phase_mean');
		sp.plot('S.rot.x.bandpass_mean');
		end
		ylim([0,1.05*sp.stat.Sc.x.clip*fc(kdx)]);
		title('Banded perpendicular','interpreter','latex');
	else
		sp.plot('S.radial.hat');
		hold on
		if (0)
		sp.plot('S.radial.brownian_phase_mean');
		sp.plot('S.radial.bandpass_mean');
		end
		ylim([0,1.15*sp.stat.Sc.radial.clip*fc(kdx)]);
		title('Spotted radial','interpreter','latex');
	end
	set(gca,'colororder',meta.colormap);
	hold on;
	xlim([0,3.5]);
	hold on
	if (0)
	legend('empirical' ...
		,['BM R^2 = ',num2str(round(sp.stat.fit.x.brownian_phase_mean.stat.r2,3))] ...
		,['BP R^2 = ',num2str(round(sp.stat.fit.x.bandpass_mean.stat.r2,3))] ...
		);
	end
	% plot density in the direction perpendicular to pattern, if pattern is banded
	splitfigure([2,4],[jdx,8],fflag,'',10.^kdx,ut,tp);
	cla();
	if (sp.stat.isisotropic)
		sp.plot('S.rot.y.hat');
		hold on
		if (0)
		sp.plot('S.rot.y.brownian_phase_across');
		end
		ylim([0,1.05*sp.stat.Sc.x.clip*fc(kdx)]);
	end
	hold on;
	xlim([0,3.5]);
	hold on
	if (0)
	legend('empirical',['BP R^2=', num2str(roundn(sp.stat.fit.y.brownian_phase_across.stat.r2,3))]);
	end
	title('Banded parallel','interpreter','latex');

	end % for kdx
	sp_a(jdx) = sp;
	end % for jdx
	if (meta.pflag)
		ps = 4;
		aspect = 4/3;
		fmt = 'pdf';

		pdfprint(11,'img/pattern-2d-band.pdf',ps,aspect,[]);
		pdfprint(12,'img/pattern-2d-band-periodogram.pdf',ps,aspect,fmt);
		pdfprint(13,'img/pattern-2d-band-autocorrelation.pdf',ps,aspect,[]);
		pdfprint(14,'img/pattern-2d-band-spectral-density.png',ps,aspect,fmt);

		aspect_ = 5/4;
		pdfprint(17,'img/pattern-2d-band-spectral-density-Sx.pdf',ps,aspect_,fmt);
		pdfprint(18,'img/pattern-2d-band-spectral-density-Sy.pdf',ps,aspect_,fmt);


	%	pdfprint(101,'img/pattern-2d-band-idealized.pdf',ps,aspect,[]);
	%	pdfprint(102,'img/pattern-2d-band-idealized-periodogram.png',ps,aspect,fmt);
	%	pdfprint(103,'img/pattern-2d-band-idealized-autocorrelation.pdf',ps,aspect,[]);
	%	pdfprint(104,'img/pattern-2d-band-idealized-spectral-density.png',ps,aspect,fmt);

	%	figure(14);
	%	colormap(gray(round(8*1.3)));
	%	clim(max(S(:))*[-0.3,1]) 
	%	pdfprint(14,'img/pattern-2d-spectral-density-grey.png',ps,aspect,fmt);
		
	%	pdfprint(15,'img/pattern-2d-transect.pdf',ps,aspect,[]);
	%	pdfprint(17,'img/pattern-2d-transect-autocorrelation.pdf',ps,aspect,[]);
	%	pdfprint(18,'img/pattern-2d-transect-spectral-density.pdf',ps,aspect,[]);

		pdfprint(21,'img/pattern-2d-spot.pdf',ps,aspect,[]);
		pdfprint(22,'img/pattern-2d-spot-periodogram.pdf',ps,aspect,fmt);
		pdfprint(24,'img/pattern-2d-spot-spectral-density.png',ps,aspect,fmt);

		pdfprint(27,'img/pattern-2d-spot-spectral-density-Sr.pdf',ps,aspect_,fmt);

	%	pdfprint(201,'img/pattern-2d-spot-idealized.pdf',ps,aspect,[]);
	%	pdfprint(202,'img/pattern-2d-spot-idealized-periodogram.png',ps,aspect,fmt);
	%	pdfprint(204,'img/pattern-2d-spot-idealized-spectral-density.png',ps,aspect,fmt);
	%	pdfprint(1000,'img/pattern-2d-band-radial-periodogram.pdf',ps);
	end
	save('mat/observed-patterns-2d.mat','sp_a');	
end % plot_observed_pattern_2d

