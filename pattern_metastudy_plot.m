%% 2022-12-02 00:25:33.449376625 +0100
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
%% plot and tabulate summary information of the meta-analysis
%
function [tab,tab_sum] = pattern_metastudy_plot(meta)
	
	if (nargin()<1)
		meta = pattern_metastudy_metadata();
	end
	fflag = meta.pflag;

	p_test  = meta.p_test;
	% note, geomean is problematic, as the resampling of S leads to zero and
	% infinitesimal small values
	% note that median and geomean require also renormalization
	avgfun  = @nanmean;

	cm = colormap_krb();
	cm = cm([3,2,1],:); % b,r,k

	% declare loaded variables
	stat =struct();
	file_C = {}; 
	Si = [];
	Ri = [];
	xi = [];
	fi = [];
	load(meta.filename.metastudy);

	% extract propertries from structures into arrays	
	Sa = cell2mat(arrayfun(@(S) rvec(S.angular_p.pdf.hp), Si,'uniformoutput',false));
	Sr = (arrayfun(@(S) rvec(S.radial.pdf.hp), Si,'uniformoutput',false));
	Sr = cell2mat(Sr);
	Sx = cell2mat(arrayfun(@(S) rvec(S.xp.pdf.hp), Si,'uniformoutput',false));
	Sy = cell2mat(arrayfun(@(S) rvec(S.y.pdf.hp), Si,'uniformoutput',false));
	Rr = cell2mat(arrayfun(@(R) rvec(R.radial), Ri,'uniformoutput',false));
%	Rx = cell2mat(arrayfun(@(R) rvec(R.x), Ri,'uniformoutput',false));
	Le_r = cvec(arrayfun(@(x) x.L_eff.r,stat));
	Le_x = cvec(arrayfun(@(x) x.L_eff.x,stat));
	Le_y = cvec(arrayfun(@(x) x.L_eff.y,stat));
	Le_aniso = sqrt(Le_x.*Le_y);

	Sxpc = cvec(cell2mat(arrayfun(@(s) rvec(s.Sc.xp.con), stat,'uniformoutput',false)));
	Syc = cvec(cell2mat(arrayfun(@(s) rvec(s.Sc.y.con), stat,'uniformoutput',false)));
	Src = cvec(cell2mat(arrayfun(@(s) rvec(s.Sc.radial.con), stat,'uniformoutput',false)));
	Stpc = cvec(cell2mat(arrayfun(@(s) rvec(s.Sc.angular_p.con), stat,'uniformoutput',false)));
	%ct = cvec(cell2mat(arrayfun(@(s) rvec(s.Sc.angular_resampled.pdf.hp), stat,'uniformoutput',false)));
	exclude      = cvec(arrayfun(@(x) x.exclude,stat));
	lxc          = cvec(arrayfun(@(x) 1./x.fc.x.hp,stat));
	lrc          = cvec(arrayfun(@(x) 1./x.fc.radial.hp,stat));
	%isisotropic  = cvec(arrayfun(@(x) x.isisotropic,stat));
	isisotropic  = cvec(arrayfun(@(x) x.isisoman,stat));
	ismodel      = cvec(arrayfun(@(x) x.ismodel,stat));
	is1d         = cvec(arrayfun(@(x) x.is1d,stat));
	isstochastic = cvec(arrayfun(@(x) x.isstochastic,stat));
	hassdf        = cvec(arrayfun(@(x) x.issdf,stat));
	p_periodic   = cvec(arrayfun(@(x) double(x.p_periodic),stat));
	fhp          = cvec(arrayfun(@(x) x.fhp,stat));

	isisotropic = logical(isisotropic);

	% relative domain size
	Le_rel_aniso = Le_aniso./lxc;
	Le_rel_r     = Le_r./lrc;
	Le_rel       = Le_rel_aniso;
	Le_rel(isisotropic) = Le_rel_r(isisotropic);

	% characteristic wavelength
	lc              = lxc;
	lc(isisotropic) = lrc(isisotropic);

	% density maxima
	S1c = Sxpc;	
	S2c = Syc;
	fdx = isisotropic==1;
	S1c(fdx) = Src(fdx);
	S2c(fdx) = Stpc(fdx).*lc(fdx);

	% regularity
	regularityx   = Sxpc./lxc;
	regularityr   = Src./lrc;
	regularity     = regularityx;
	regularity(isisotropic) = regularityr(isisotropic);
	fc = 1./lc;

	% check sanity of extimates
	% the analysis fails for pattern without a characteristic length scale,
	% these have to be excluded
	fdx = (lc == 0) | ~isfinite(lc);
	exclude(fdx) = 1;
	fdx = find(fdx & ~exclude);
	printf('Number of non-exluded patterns where lc = 0: %d\n',length(fdx));
	disp(fdx)

	fdx = ~isfinite(p_periodic);
	exclude(fdx) = 1;
	fdx = find(fdx & ~exclude);
	printf('Number of non-exluded patterns where p_periodic is not finite: %d\n',length(fdx));
	disp(fdx)

	fdx = (fc<fhp & ~exclude);
	exclude(fdx) = 1;
	fdx = find(fdx & ~exclude);
	printf('Number of non-exluded patterns where fc < fhpe: %d\n',length(fdx));

	fdx = isnan(ismodel);
	exclude(fdx) = 1;
	printf('Number of patterns where model/nature was not specified:',find(fdx));
	ismodel(isnan(ismodel)) = false;

	
	% test for differences between natural and model generated patterns
	fdx0 =  (hassdf) &  ismodel & (~is1d) & (~isstochastic) & ~exclude;
	fdx1 =  (hassdf) & ~ismodel & (~is1d) & (~isstochastic) & ~exclude;
	if (exist('mediantest','file'))
		p = mediantest(regularity(fdx0),regularity(fdx1));
		printf('Test model vs nature for different median %e\n',p);
	else
		printf('mediantest not installed, skipping test\n');
	end
	
	% generate tex file for table of analysis results
	fid = fopen('mat/metastudy-table.tex','w');

	% generate table with summary statistics
	nL  = 10;
	lc_ = inf;
	tab_sum = table;
	fdx = (exclude == 0);
	fprintf('N total %d\n',sum(fdx));
	id = 1;
	tab_sum{id,1} = "Total";
	tab_sum{id,2} = sum(fdx);
	tab_sum{id,3:5} = NaN;

	% confidence interval
	%pc = 0.64;
	pc = 1-2*(1-normcdf(1))	
	%pc = 1-2*(1-normcdf(2))	

	fdx = (hassdf == 1) & (~exclude) & (nL*lc < lc_) & (ismodel == 0);
	printf('Nature all: %d %0.2f %0.2f %0.2f\n',sum(fdx),mean(p_periodic(fdx)<p_test), median(S1c(fdx)./lc(fdx)),median(S2c(fdx)./lc(fdx)));
	tabulate(2,'Nature all');
%	tab_sum{id,6}  = hassdf(fdx);

	fdx = (hassdf == 1) & ~exclude & (nL*lc < lc_) & (ismodel == 0) & (isisotropic == 0);
	printf('Nature aniso: %d %0.2f %0.2f %0.2f\n',sum(fdx),mean(p_periodic(fdx)<p_test), median(Sxpc(fdx)./lc(fdx)),median(Syc(fdx)./lc(fdx)));
	tabulate(3,'Nature anisotropic');

	fdx = (hassdf == 1) & (~exclude) & (nL*lc < lc_) & (ismodel == 0) & (isisotropic == 1);
	printf('Nature   iso: %d %0.2f %0.2f %0.2f\n',sum(fdx), mean(p_periodic(fdx)<p_test), median(Src(fdx)./lc(fdx)),median(Stpc(fdx)));
	tabulate(4,'Nature isotropic');

	fdx = (hassdf == 1) & (~exclude) & (nL*lc < lc_) & (cvec(isstochastic) == 0) & (ismodel == 1) & (is1d == 0); 
	fprintf('2D Model Total: %d %0.2f %0.2f %0.2f\n',sum(fdx),mean(p_periodic(fdx)<p_test), median(S1c(fdx)./lc(fdx)),median(S2c(fdx)./lc(fdx)));
	tabulate(5,'2D model deterministic all');

	fdx = (hassdf == 1) & (exclude == 0) & (nL*lc < lc_) & (cvec(isstochastic) == 0) & (ismodel == 1) & (isisotropic == 0) & (is1d == 0); 
	fprintf('2D Model Aniso: %d %0.2f %0.2f %0.2f\n',sum(fdx),mean(p_periodic(fdx)<p_test), median(Sxpc(fdx)./lc(fdx)),median(Syc(fdx)./lc(fdx)));
	tabulate(6,'2D model deterministic anisotropic');

	fdx = (hassdf == 1) & (exclude == 0) & (nL*lc < lc_) & (cvec(isstochastic) == 0) & (ismodel == 1) & (isisotropic == 1) & (is1d == 0);
	fprintf('2D model iso : %d %0.2f %0.2f %0.2f\n',sum(fdx), mean(p_periodic(fdx)<p_test), median(Src(fdx)./lc(fdx)),median(Stpc(fdx)));
	tabulate(7,'2D model deterministic isotropic');

	fdx = (hassdf == 1) & (exclude == 0) & (nL*lc < lc_) & (cvec(isstochastic) == 1) & (ismodel == 1) & (isisotropic == 1) & (is1d == 0);
	fprintf('2D Model stoch iso: %d %0.2f %0.2f %0.2f\n',sum(fdx), mean(p_periodic(fdx)<p_test),median(Src(fdx)./lc(fdx)),median(Stpc(fdx)));
	tabulate(8,'2D model stochastic isotropic');

	fdx = (hassdf == 0) & (exclude == 0) & (nL*lc < lc_) & (cvec(isstochastic) == 0) & (ismodel == 1) & (isisotropic == 1) & (is1d == 0);
	fprintf('2D Model, w/o sdf: %d %0.2f %0.2f %0.2f\n',sum(fdx), mean(p_periodic(fdx)<p_test),median(Src(fdx)./lc(fdx)),median(Stpc(fdx)));
	tabulate(9,'2D model w/o sdf');

	fdx = (hassdf == 1) & (exclude == 0) & (nL*lc < lc_) & (cvec(isstochastic) == 0) & (ismodel == 1) & (is1d == 1);
	fprintf('1D Model: %d %0.2f %0.2f\n',sum(fdx),mean(p_periodic(fdx)<p_test),median(Sxpc(fdx)./lc(fdx)));
	tabulate(10,'1D model deterministic');
	disp(tab_sum)
	tab_sum.Properties.VariableNames = {'Group','N','$p<$0.05','$S_{1c}/\lambda_c$','l1','u1','$S_{2c}/\lambda_c$','l2','u2','$L_\{eff}$'};
	tab_sum.Properties.RowNames = tab_sum.Group;
	disp(tab_sum);

	fid2 = fopen('mat/metastudy-table-summary.tex','w');
	s = table2tex(tab_sum,2);
	fprintf(fid2,s);
	fclose(fid2);

	% detailed table pattern by pattern
	tab = table();

	j = 0;
	for idx=1:length(regularity)
	% tex-entry
		f = file_C{idx};
		d = dirname(f);
		b = basename(f);
		figid = regexprep(b,'_.*','');
		L_eff = NaN;
		f = f(1:end-4);
		f = regexprep(f,'patterns/metastudy','metastudy');
		f = ['',dirname(f),'/crop/',basename(f)];
		%suffix1 = {'_{x}^+','_r'};
		%suffix2 = {'_y','_{\theta}^+'};
		%subsrcipt1 = 'xr';
		%superscript1 = '+ ';
		%subsrcipt2 = {'y','\theta'};
		%superscript2 = ' 
		root = 'img/metastudy-result/';
		if (~exclude(idx))
		j = j+1;
		if (~is1d(idx))
			img2d = sprintf('	  \\includegraphics[height=0.22\\textwidth]{%s/%s-density-2d-4-crop.pdf}',root,f); ...
			imgy  = sprintf('	  \\includegraphics[height=0.22\\textwidth]{%s/%s-density-Sy-4-crop.pdf}',root,f); ...
			if (isisotropic(idx))
				density_1_label = 'S_r';
				density_2_label = 'S_{\theta}^+';
				reg_1_label = 'S_{rc}/\lambda_c';
				reg_2_label = 'S_{{\theta}c}^+';
				reg_2 = stat(idx).Sc.angular_p.con;
			else
				density_1_label = 'S_x^+';
				density_2_label = 'S_y';
				reg_1_label = 'S_{xc}^+/\lambda_c';
				reg_2_label = 'S_y/\lambda_c';
				reg_2 = stat(idx).Sc.y.con/lc(idx);
			end
		else
			img2d = 'X (1d)';
			imgy  = 'X (1d)';
			reg_1_label = 'S_{xc}^+/\lambda_c';
			reg_2_label = 'N/A';
			reg_2       = NaN;
		end
		fprintf(fid,[...
		'\\begin{tabular}{cccc|c|c|c}\n' ...
		'\\multicolumn{6}{l}{\\parbox{\\textwidth}{%s %s %s}}\n' ...
		'\\\\\\hline\n' ...
		'Density $S_{xy}$ & Density $%s$ & Density $%s$ & $%s$ & $%s$ & $L_{eff}/\\lambda_c$ & $p$ \n' ...
		],...	
		num2str(idx),...
		figid, ...
		basename(d),...
		density_1_label, ...
		density_2_label, ...
		reg_1_label, ...
		reg_2_label ...
		);
		fprintf(fid,[...	
		   '\\\\\\hline\n' ...
		'  %s\n' ...
		'& \\includegraphics[height=0.22\\textwidth]{%s/%s-density-Sx-4-crop.pdf}\n' ...
		'& %s\n'], ...
		img2d,root,f,imgy);
	
		fprintf(fid,[ ...
		'& \\raisebox{11ex}{%0.2f}\n' ...
		'& \\raisebox{11ex}{%0.2f}\n' ...
		'& \\raisebox{11ex}{%0.1f}\n' ...
		'& \\raisebox{11ex}{%0.3f}\n' ...
		... '\\end{tblr}' 
		'\\end{tabular}\n\\\\' ...
		],...
		regularity(idx),reg_2,Le_rel(idx),p_periodic(idx));
	
		% csv-entry
		dirname_  = dirname(file_C{j});
		word_C    = strsplit(dirname_,'-');
		try
			tab.author{j} = word_C{end-1};
			tab.year(j)   = str2double(word_C{end});
		catch e
			tab.author{j} = 'this study';
			tab.year(j) = NaN;
		end
		switch (isisotropic(j))
		case {1}
			tab.Isotropy{j} = 'isotropic';
		case {0}
			tab.Isotropy{j} = 'anisotropic';
		otherwise
			tab.Isotropy{j} = '/';
		end % switch
		tab.regularity_x(j) = regularity(j);
		tab.regularity_y(j) = reg_2;
		tab.Le_rel(j)       = Le_rel(j);
		tab.p_periodic(j)   = p_periodic(j);
		switch (is1d(j))
		case {1}
			tab.dimension(j) = 1;
		case {0}
			tab.dimension(j) = 2;
		otherwise
			tab.dimension(j) = NaN;
		end % switch
		if (ismodel(j)==1)
			tab.type{j}      = 'model';
		if (isstochastic==1)
			tab.mtype{j} = 'heterogeneous';
		else
			tab.mtype{j} = 'homogeneous';
		end % else of isstoch
		else
			tab.type{j} = 'nature';
			tab.mtype{j} = '';
		end % else of ismodel
		tab.hassdf{j} = hassdf(j);

	end % if ~exclude
	
	end % for idx
	fclose(fid);
	filename = meta.filename.patterns_literature_stat_csv;
	writetable(tab,filename);
	copyfile(filename,[filename(1:end-4),'-',datestr(now(),'yyyy-mm-dd'),'.csv']);

	for xy=0:1
	for isiso_=0:1
	% plot density	
	splitfigure([2,2],[30,1+isiso_+2*xy],fflag);
	cla();
	
	for ismodel_=0:1
	if (ismodel_)
		fdx =  (hassdf) &  ismodel & (~is1d) & (isisotropic == isiso_) & (~isstochastic) & ~exclude;
	else
		fdx =  (hassdf) & ~ismodel & (isisotropic == isiso_) & ~exclude;
	end % if jdx
	
	if (xy==0)
	if (ismodel_ == 0)
		fdx0 = fdx;
	else
		p = mediantest(regularity(fdx0),regularity(fdx));
		fprintf('Test iso = %d model vs nature for different median %e\n',isiso_,p);
	end
	end

	if (isiso_)
		% note that the resampled Sc is already normalized by lc
		if (xy)
			S  = cvec(avgfun(Sa(fdx,:)));
			f  = cvec(fi.angular);
		else
			S = cvec(avgfun(Sr(fdx,:)));
			f = cvec(fi.radial);
		end
		plot(f,S,'linewidth',1);
		if (xy)
			xlim([-1,1]*pi/2);
			xlabel('Angle $\theta$','interpreter','latex');
			ylabel('Density $S_\theta^+$','interpreter','latex');
			set(gca,'xtick',[-1/2,-1/4,0,1/4,1/2]*pi,'xticklabel',{'-\pi/2','-\pi/4','0','\pi/4','\pi/2'});
		else
			xlim([0,2.5]);
			xlabel('Wavenumber $k_r / k_c$','interpreter','latex');
			ylabel('Density $S_r/\lambda_c$','interpreter','latex');
		end
	else
		if (xy)
			S = avgfun(Sy(fdx,:));
			f = fi.y;
		else
			S = avgfun(Sx(fdx,:));
			f = fi.x;
		end
		plot(f,S,'linewidth',1);
		if (xy)
		xlabel('Wavenumber $k_y/k_c$','interpreter','latex');
		ylabel('Density $S_y / \lambda_c$','interpreter','latex');
		xlim([-2,2]);
		else
			xlim([0,2.5]);
			xlabel('Wavenumber $k_x / k_c$','interpreter','latex');
			ylabel('Density $S_x^+/\lambda_c$','interpreter','latex');
			legend('Nature','Model');
		end
	end % else of isiso
		hold on
	end % for ismodel_
		%colormap('default');
		%colormap(cm([4,2],:));
		set(gca,'colororder',cm);
	end % for xy
	
	end % for isiso
	
	% coorelation analysis
	fdx = (hassdf) & (ismodel==0) & (~is1d) & (isisotropic == 0) & (~isstochastic) & ~exclude;
	c(1) = kendall_to_pearson(corr(cvec(Sxpc(fdx)),cvec(Syc(fdx))));

	splitfigure([2,2],[100,3],fflag);
	cla();
	plot(log10(cvec(Sxpc(fdx))./lc(fdx)),log10(cvec(Syc(fdx)./lc(fdx))),'.')
	hold on
	fdx = hassdf & ismodel & (~is1d) & (isisotropic == 0) & (~isstochastic) & ~exclude;
	c(2) = kendall_to_pearson(corr(cvec(Sxpc(fdx)),cvec(Syc(fdx))));
	plot(log10(cvec(Sxpc(fdx))./lc(fdx)),log10(cvec(Syc(fdx)./lc(fdx))),'.')
	fdx = hassdf & (ismodel==0) & (~is1d) & (isisotropic == 1) & (~isstochastic) & ~exclude;
	c(3) = kendall_to_pearson(corr(cvec(Src(fdx))./lc(fdx),cvec(Stpc(fdx))));

	splitfigure([2,2],[100,4],fflag);
	cla();
	plot(log10(cvec(Src(fdx))./lc(fdx)),cvec(Stpc(fdx)),'.')
	hold on
	fdx = hassdf & ismodel & (~is1d) & (isisotropic == 1) & (~isstochastic) & ~exclude;
	c(4) = kendall_to_pearson(corr(cvec(Src(fdx))./lc(fdx),cvec(Stpc(fdx))));
	plot(log10(cvec(Src(fdx))./lc(fdx)),cvec(Stpc(fdx)),'.')

	printf('Correlation nature aniso Sxpc,Syc: %f\n',c(1));
	printf('Correlation model  aniso Sxpc,Syc: %f\n',c(2));
	printf('Correlation nature   iso Src/lc,Stpc: %f\n',c(3));
	printf('Correlation model    iso Src/lc,Stpc: %f\n',c(4));

	% plot correlation
	splitfigure([2,2],[100,1],fflag);
	cla();
	for idx=1:2
		plot(idx,c(idx),'*','color',cm(idx,:))
		hold on
	end
	ax = gca;
	set(ax(1),'ylim',[-0.3,1.05])
	ylabel(ax(1),'corr($S_{xc}^+$,$S_{yc}$)','interpreter','latex');
	xlim([0.5,2.5]);
	set(gca,'xtick',1:4,'xticklabel',{'nature','model','nature','model'},'xticklabelrot',45);
	grid on
	daspect([2,1,1]);
%	yyaxis right


	% plot correlation
	splitfigure([2,2],[100,2],fflag);
	for idx=3:4
		plot(idx-2,c(idx),'*','color',cm(idx-2,:))
		hold on
	end
	ax=gca
	%set(ax,'ycolor','k');
	%set(ax(1),'ylim',[-0.2,1.09])
	set(ax(1),'ylim',[-0.3,1.05])
	ylabel(ax,'corr($S_{rc}$,$S_{sc}$)','interpreter','latex');
	xlim([0.5,2.5]);
%	text(0.6,0.075,'Anisotropic');
%	text(2.8,0.075,'Isotropic');
%	vline(2.5,'linestyle','--','color','k');
	set(gca,'xtick',1:4,'xticklabel',{'nature','model','nature','model'},'xticklabelrot',45);
	grid on
	daspect([2,1,1])	

	% plot autocorrelation function	
	for isiso_=0:1
	splitfigure([2,2],[40,1+isiso_],fflag);
	cla();
	for ismodel_=0:1
	if (ismodel_)
		fdx = ismodel & (~is1d) & (isisotropic == isiso_) & (~isstochastic) & ~exclude;
	else
		fdx = ~ismodel & (isisotropic == isiso_) & ~exclude;
	end % if jdx
		
	if (isiso_)
if (0)
		R_ = mean(Rr(fdx,:));
		plot(xi,R_,'linewidth',1);
		hold on
		xlim([0,2.5]);
		ylabel('R_r');
		xlabel('Lag distance r/\lambda_c');
		if (ismodel_ && isiso_)
		o = 2*pi*xi;
	%	plot(xi,besselj(0,2*pi*xi));
		plot(xi,sqrt(2/pi).*1./sqrt(o),'k--');
		h=plot(xi,-sqrt(2/pi).*1./sqrt(o),'k--');
		h.HandleVisibility='off';
		legend('nature','model','$\hat J_0$','interpreter','latex');
		ylim([-1.05,1.05]);
		end
end
	else
if (0)
		%R_ = mean(Rx(fdx,:));
		%plot(xi,R_,'linewidth',1);
		hold on
		xlim([0,2.5]);
		ylabel('R_x');
		xlabel('Lag distance x/\lambda_c');
		if (ismodel_)
		plot(xi,ones(size(xi)),'k--');
		h=plot(xi,-ones(size(xi)),'k--');
		end
		ylim([-1.05,1.05]);
end
	end
		grid on
		
	end % ismodel_
	end % isiso
	
	% plot density

	isisotropic_        = [0, 0, 1, 1, 1,   0];
	is1d_               = [0, 0, 0, 0, 0,   1];
	ismodel_            = [0, 1, 0, 1, 1,   1];
	isstochastic_       = [0, 0, 0, 0, 1,   0];
	
	leg_C = {'nature','model','nature','deterministic','stochastic','model'};
	
	id = 1:length(isisotropic_);
	
	q=[];
	np = [];
	for idx=id
		fdx = (   (~exclude) ...
		        & (ismodel  == ismodel_(idx)) ...
		        & (is1d     == is1d_(idx) | ~ismodel_(idx)) ...
		        & ((isstochastic  == isstochastic_(idx)) | ~ismodel_(idx)) ...
		        & (isisotropic   == isisotropic_(idx) | is1d_(idx) == 1) ...
		      );
		np(idx,1)       = sum(fdx);
		q(:,idx)        = quantile(regularity(fdx),[0.25,0.5,0.75]);
		qL_rel(:,idx)   = quantile(Le_rel(fdx),[0.25,0.5,0.75]);
	end
	da=[1.5   15.875    1];
	da=[1.5   10    1];
	da_=[1.5   0.375    1];

	% plot quantiles of the regularity
	splitfigure([2,2],[10,1],fflag);
	cla();
	id_ = [1:2,6];
	idp = [1:3];
	h=errorbar(1,q(2,1),q(2,1)-q(1,1),q(3,1)-q(2,1),'*','color',cm(1,:));
	h.HandleVisibility='off';
	hold on
	h=errorbar(2,q(2,2),q(2,2)-q(1,2),q(3,2)-q(2,2),'*','color',cm(2,:));
	h.HandleVisibility='off';
	h=errorbar(3,q(2,6),q(2,6)-q(1,6),q(3,6)-q(2,6),'k*');
	h.HandleVisibility='off';
	% disp(sum(np))
	k = 0;
	for idx=id_
		k=k+1;
		text(k,q(2,idx),[' ',num2str(np(idx))]);
	end
	xlim([min(idp)-0.5,max(idp)+0.5]);
	set(gca,'yscale','log');
	ylabel('Regularity $S_{xc}^+/\lambda_c$','interpreter','latex');
	set(gca,'xtick',idp,'xticklabel',{'nature','2D-model','1D-model'},'xticklabelrot',45);
	set(gca,'ytick',2.^(-2:6));
	ylim([0.4,20]) 
	vline(2.5,'linestyle','--','color','k');
	ty = 1.4*32-20;
	%text(1-0.67*0-0.25,ty,'2D')
	%text(3-0.0,ty,'1D')
	daspect(da)

	% plot ranges of the primary regularity with ci
		splitfigure([2,2],[50,1],fflag);
		cla();
		r = {'Nature anisotropic','2D model deterministic anisotropic','1D model deterministic'}
		c = {'l1','$S_{1c}/\lambda_c$','u1'}
		for idx=1:length(r)
		errorbar(idx,tab_sum{r{idx},c{2}},tab_sum{r{idx},c{2}}-tab_sum{r{idx},c{1}},tab_sum{r{idx},c{3}}-tab_sum{r{idx},c{2}},'*','color',cm(idx,:));
		hold on
		text(idx,tab_sum{r{idx},c{2}},[' ',num2str(tab_sum{r{idx},'N'})]);
		end
		ylim([0.4,20]);
		xlim([0.5,3.5]);
		set(gca,'yscale','log');
		ylabel('Regularity $S_{xc}^+/\lambda_c$','interpreter','latex');
	set(gca,'ytick',2.^(-2:6));
		id_ = 1:3;
		set(gca,'xtick',id_,'xticklabel',{'nature','2D-model','1D-model'},'xticklabelrot',45);
	 	vline(2.5,'linestyle','--','color','k');
	daspect(da);

		splitfigure([2,2],[50,2],fflag);
		cla
		r = {'Nature isotropic','2D model deterministic isotropic','2D model stochastic isotropic'}
		c = {'l1','$S_{1c}/\lambda_c$','u1'}
		for idx=1:length(r)
		errorbar(idx,tab_sum{r{idx},c{2}},tab_sum{r{idx},c{2}}-tab_sum{r{idx},c{1}},tab_sum{r{idx},c{3}}-tab_sum{r{idx},c{2}},'*','color',cm(idx,:));
		hold on
		text(idx,tab_sum{r{idx},c{2}},[' ',num2str(tab_sum{r{idx},'N'})]);
		end
		ylim([0.4,20]);
		xlim([0.5,3.5]);
		set(gca,'yscale','log');
		ylabel('Regularity $S_{rc}/\lambda_c$','interpreter','latex');
	set(gca,'ytick',2.^(-2:6));
		id_ = [1,2,2.3,3,3.3];
		set(gca,'xtick',id_,'xticklabel',{'nature','determinis-','tic model','stochastic','model'},'xticklabelrot',45);
	daspect(da);

		splitfigure([2,2],[50,3],fflag);
		cla();
		r = {'Nature anisotropic','2D model deterministic anisotropic'} %,'1D model deterministic'}
		c = {'l2','$S_{2c}/\lambda_c$','u2'}
		for idx=1:length(r)
		errorbar(idx,tab_sum{r{idx},c{2}},tab_sum{r{idx},c{2}}-tab_sum{r{idx},c{1}},tab_sum{r{idx},c{3}}-tab_sum{r{idx},c{2}},'*','color',cm(idx,:));
		hold on
		text(idx,tab_sum{r{idx},c{2}},[' ',num2str(tab_sum{r{idx},'N'})]);
		end
		ylim([0.4,20]);
		xlim([0.5,2.5]);
		set(gca,'yscale','log');
		ylabel('Regularity $S_{yc}/\lambda_c$','interpreter','latex');
	set(gca,'ytick',2.^(-2:6));
		id_ = 1:3;
		set(gca,'xtick',id_,'xticklabel',{'nature','2D-model'},'xticklabelrot',45);
	daspect(da);


		splitfigure([2,2],[50,4],fflag);
		cla
		r = {'Nature isotropic','2D model deterministic isotropic','2D model stochastic isotropic'}
		c = {'l2','$S_{2c}/\lambda_c$','u2'}
		for idx=1:length(r)
		errorbar(idx,tab_sum{r{idx},c{2}},tab_sum{r{idx},c{2}}-tab_sum{r{idx},c{1}},tab_sum{r{idx},c{3}}-tab_sum{r{idx},c{2}},'*','color',cm(idx,:));
		hold on
		text(idx,tab_sum{r{idx},c{2}},[' ',num2str(tab_sum{r{idx},'N'})]);
		end
		ylim([0.25,1]);
		xlim([0.5,3.5]);
		set(gca,'yscale','log');
		ylabel('Regularity $S_{{\theta}c}^+/\lambda_c$','interpreter','latex');
	set(gca,'ytick',2.^(-2:6));
		id_ = [1,2,2.3,3,3.3];
		set(gca,'xtick',id_,'xticklabel',{'nature','determinis-','tic model','stochastic','model'},'xticklabelrot',45);
	daspect(da_);
	%end 

	% plot regularity with interquartile ranges
	% isotropic
	splitfigure([2,2],[10,2],fflag);
	cla();
	id_=3:5;
	h=errorbar(3,q(2,3),q(2,3)-q(1,3),q(3,3)-q(2,3),'*','color',cm(1,:));
	hold on
	h=errorbar(4,q(2,4),q(2,4)-q(1,4),q(3,4)-q(2,4),'*','color',cm(2,:));
	h=errorbar(5,q(2,5),q(2,5)-q(1,5),q(3,5)-q(2,5),'k*');
	%vline(4.5,'linestyle','--','color','k');
	for idx=id_
		text(idx,q(2,idx),[' ',num2str(np(idx))]);
	end
	xlim([min(id_)-0.5,max(id_)+0.5]);
	set(gca,'yscale','log');
	ylabel('Regularity $S_{rc}/\lambda_c$','interpreter','latex');
	set(gca,'xtick',id,'xticklabel',leg_C,'xticklabelrot',45);
	set(gca,'ytick',2.^(-2:6));
	ylim([0.4,20]) 
	id_ = [3,4,4.3,5,5.3];
	set(gca,'xtick',id_,'xticklabel',{'nature','determinis-','tic model','stochastic','model'},'xticklabelrot',45);
	daspect(da);

	% plat regularity (combine isotropic and anisotropic models)	
	splitfigure([2,2],[10,3],fflag);
	cla
	q   = [];
	sto = [0,0,1];
	mo  = [0,1,1];
	for idx=1:3
		lc=cvec(lc);
		fdx       = ~isisotropic & ~exclude & (ismodel == mo(idx)) & (isstochastic==sto(idx)) & (is1d == 0);
		q(:,idx)  = quantile(cvec(Syc(fdx))./lc(fdx),[0.25,0.5,0.75]);
		fdx       = isisotropic & ~exclude & (ismodel == mo(idx)) & (isstochastic==sto(idx)) & (is1d == 0);
		qt(:,idx) = quantile(cvec(Stpc(fdx)),[0.25,0.5,0.75]);
	end
	for idx=1:2
		h=errorbar(idx,q(2,idx),q(2,idx)-q(1,idx),q(3,idx)-q(2,idx),'*','color',cm(idx,:));
		hold on
	end
	%h=errorbar(3,q(2,3),q(2,3)-q(1,3),q(3,3)-q(2,3),'k*');
	xlim([0.5,2.5]);
	set(gca,'yscale','log');
	ylabel('Regularity $S_{yc}/\lambda_c$','interpreter','latex');
	set(gca,'xtick',1:2,'xticklabel',{'nature','model'},'xticklabelrot',45);
	ylim([0.4,20]) 
	set(gca,'ytick',[2.^(-2:1:5)])
	daspect(da)
	
	% angular
	splitfigure([2,2],[10,4],fflag);
	cla
	for idx=1:3
		h=errorbar(idx,qt(2,idx),qt(2,idx)-qt(1,idx),qt(3,idx)-qt(2,idx),'*','color',cm(idx,:));
		hold on
	end
	xlim([0.5,3.5]);
	set(gca,'yscale','log');
	ylabel('Regularity $S_{{\theta}c}^+$','interpreter','latex');
	% = \frac{S_{cs}}{\lambda_c}$','interpreter','latex');
	%set(gca,'xtick',1:3,'xticklabel',{'nature','deterministic','stochastic'},'xticklabelrot',45);
	id_ = [1,2,2.3,3,3.3];
	set(gca,'xtick',id_,'xticklabel',{'nature','determinis-','tic model','stochastic','model'},'xticklabelrot',45);
	%ylim([0.4,20]) 
	%d_=[1.5   0.375    1];
	%d_=[1.5   4.6*0.375    1];
	d_=[1.5  0.375    1];
	%ylim([1,4.0]) 
	ylim([0.25,1]) 
	set(gca,'ytick',[2.^(-2:1:5)])
	daspect(da_)
	%daspect(d)

	% plot regularity, simple	
	ismodel_       = [0, 1];
	isstochastic_  = [0, 0];
	splitfigure([2,2],[20,1],fflag);
	cla();
	q=[];
	np = [];
	id = 1:length(ismodel_);
	for idx=id
		fdx = (   (0 == exclude) ...
		        & (ismodel  == ismodel_(idx)) ...
		        & (cvec(isstochastic)  == isstochastic_(idx)) ....
			& ~exclude ...
		      );
		np(idx,1) = sum(fdx);
		q(:,idx) = quantile(regularity(fdx),[0.25,0.5,0.75]);
	end
	 h=errorbar(id,q(2,:),q(2,:)-q(1,:),q(3,:)-q(2,:),'*');
	disp(sum(np))
	for idx=id
		text(idx,q(2,idx),[' ',num2str(np(idx))]);
	end
	 xlim([min(id)-0.5,max(id)+0.5]);
	 h.HandleVisibility='off';
	 set(gca,'yscale','log');
	 ylabel('Regularity S_c/\lambda_c');
	 set(gca,'xtick',id,'xticklabel',leg_C,'xticklabelrot',45);
	 set(gca,'ytick',2.^(-2:6));
	 ylim([0.25,6]) 
	 vline(2.5,'linestyle','--','color','k');
	 vline(5.5,'linestyle','--','color','k');
	 ty = 1.4*64;
	 text(1-0.67*0-0.25,ty,'2D')
	 text(6,ty,'1D')

	% R-simple	
	splitfigure([2,2],[60,1],fflag);
	cla();
	ismodel_ = [0,1,1];
	is1d_    = [0,0,1];
	f1d = fi.x;
	for idx=1:3
		fdx =   (ismodel_(idx) == ismodel) ...
		      & (is1d_(idx) == is1d) ...
		      & (exclude == 0) ...
	              & ~exclude;
		S = Sr;
		S(~isisotropic,:) = Sx(~isisotropic,:);
if (0)
		R = Rr;
		%R(~isisotropic,:) = Rx(~isisotropic,:);
		R1d_(:,idx) = nanmean(R(fdx,:));
end	
		S1d_(:,idx) = nangeomean(S(fdx,:));
	end
	S1d_ = S1d_./(sum(S1d_)*(f1d(2)-f1d(1)));
	
	plot(f1d,S1d_,'linewidth',1);
	xlabel('Wavenumber k/k_c');
	ylabel('Density S_c/\lambda_c');
	legend('real','model 2d','model 1d');
	xlim([0,2.5])
	set(gca,'colororder',colormap_krb());
	
	if (0)
		S1d_(~isfinite(S1d_))=0;
		m=0;
		R1d = real(ifft([S1d_; zeros(m,size(S1d_,2)); flipud(S1d_(2:end,:))]));
		R1d = R1d./R1d(1,:);
		df=f1d(2)-f1d(1);
		%L = 1./df;
		x=linspace(0,L,length(R1d))';
	end
	
	% plot autocorrelation
	splitfigure([2,2],[60,2],fflag);
	cla
if (0)
	plot(xi,R1d_./R1d_(1,:),'linewidth',1);
end
	xlabel('Lag Distance x/\lambda_c');
	ylabel('Autocorrelation R');
	xlim([0,2.25])
	set(gca,'colororder',colormap_krb());
	
	splitfigure([2,2],[20,3],fflag);
	cla();
	ismodel_ = [0,1];
	is1d_    = [0,0];
	S1d_ = [];
	for idx=1:length(ismodel_)
		fdx =   (ismodel_(idx) == ismodel) ...
		      & (exclude == 0) ...
	              & ~exclude;
	
		S = Sr;
		S(~isisotropic,:) = Sx(~isisotropic,:);
		S1d_(:,idx) = nangeomean(S(fdx,:));
if (0)
		R = Rr;
		%R(~isisotropic,:) = Rx(~isisotropic,:);
		R1d_(:,idx) = nanmean(R(fdx,:));
end	
	end
	S1d_ = S1d_./(sum(S1d_)*(f1d(2)-f1d(1)));
	
	plot(f1d,S1d_,'linewidth',1);
	xlabel('Wavenumber k/k_c');
	ylabel('Density S_c/\lambda_c');
	legend('real','model','model 1d');
	xlim([0,2.5])
	set(gca,'colororder',colormap_krb());
	
	fdx = ~exclude & ~exclude;
	f = file_C(fdx);
	f = cellfun(@dirname,f,'uniformoutput',false)';
	f =unique(f)


	if (meta.pflag)
		ps = 3.5;
		ps_ = 4;
		% plot with quartiles
		pdfprint(101,'img/metastudy-regularity-Sxc-q.pdf',ps_);
		pdfprint(102,'img/metastudy-regularity-Src-q.pdf',ps_);
		pdfprint(103,'img/metastudy-regularity-Syc-q.pdf',ps_);
		pdfprint(104,'img/metastudy-regularity-Stc-q.pdf',ps_);

		pdfprint(1001,'img/metastudy-regularity-correlation-Sxpc-Syc.pdf',ps_);
		pdfprint(1002,'img/metastudy-regularity-correlation-Src-Stc.pdf',ps_);

		% plot with confidence intervals
		pdfprint(501,'img/metastudy-regularity-Sxc-ci.pdf',ps_);
		pdfprint(502,'img/metastudy-regularity-Src-ci.pdf',ps_);
		pdfprint(503,'img/metastudy-regularity-Syc-ci.pdf',ps_);
		pdfprint(504,'img/metastudy-regularity-Stc-ci.pdf',ps_);
		
		figure(103);
		axis square
		pdfprint(103,'img/density-literature.pdf',ps);
		pdfprint(104,'img/autocorrelation-literature.pdf',ps);
	
		pdfprint(201,'img/regularity-literature-simple.pdf',ps);
		pdfprint(203,'img/density-literature-simple.pdf',ps);
	%	pdfprint(204,'img/autocorrelation-literature-simple.pdf',ps);
	
		pdfprint(301,'img/metastudy-density-Sx.pdf',ps);
		pdfprint(302,'img/metastudy-density-S-radial.pdf',ps);
		pdfprint(303,'img/metastudy-density-Sy.pdf',ps);
		pdfprint(304,'img/metastudy-density-S-angular.pdf',ps);
	
		pdfprint(401,'img/metastudy-autocorrelation-Rx.pdf',ps);
		pdfprint(402,'img/metastudy-autocorrelation-Rr.pdf',ps);
	
	
	end % if pflag

	function tabulate(id,title_)
		tab_sum{id,1} = {title_};
		tab_sum{id,2} = sum(fdx);
		tab_sum{id,3} = round(mean(p_periodic(fdx)<p_test),2);
		tab_sum{id,4} = round(median(S1c(fdx)./lc(fdx)),2);
		[me,sme,lme,rme] = median_man(S1c(fdx)./lc(fdx),pc);
		tab_sum{id,5} = lme;
		tab_sum{id,6} = rme;
		tab_sum{id,7} = round(median(S2c(fdx)./lc(fdx)),2);
		[me,sme,lme,rme] = median_man(S2c(fdx)./lc(fdx),pc);
		tab_sum{id,8} = lme;
		tab_sum{id,9} = rme;
		tab_sum{id,10} = median(Le_rel(fdx));
	end % tabulate
end % pattern_metastudy_plot


