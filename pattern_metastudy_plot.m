% 2022-12-02 00:25:33.449376625 +0100
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
% plot and tabulate summary information of the meta-analysis
%
function pattern_metastudy_plot(meta)
	
	if (nargin()<1)
		meta = pattern_analysis_metadata();
	end
	fflag = meta.pflag;

	p_test  = meta.p_test;
	avgfun  = @nangeomean;

	cm = meta.colormap;
	cm = cm([3,2,1],:); % b,r,k

	% 
	load(meta.filename.metatstudy);

	% extract propertries from structures into arrays	
	Sa = cell2mat(arrayfun(@(S) rvec(S.angular.pdf.clip), Si,'uniformoutput',false));
	Sr = (arrayfun(@(S) rvec(S.radial.pdf.clip), Si,'uniformoutput',false));
	Sr = cell2mat(Sr);
	Sx = cell2mat(arrayfun(@(S) rvec(S.x.pdf.clip), Si,'uniformoutput',false));
	Sy = cell2mat(arrayfun(@(S) rvec(S.y.pdf.clip), Si,'uniformoutput',false));
	Rr = cell2mat(arrayfun(@(R) rvec(R.radial), Ri,'uniformoutput',false));
%	Rx = cell2mat(arrayfun(@(R) rvec(R.x), Ri,'uniformoutput',false));
	Le_r = cvec(arrayfun(@(x) x.L_eff.r,stat));
	Le_x = cvec(arrayfun(@(x) x.L_eff.x,stat));

	Scx = cvec(cell2mat(arrayfun(@(s) rvec(s.Sc.x.clip), stat,'uniformoutput',false)));
	Scy = cvec(cell2mat(arrayfun(@(s) rvec(s.Sc.y.clip), stat,'uniformoutput',false)));
	Scr = cvec(cell2mat(arrayfun(@(s) rvec(s.Sc.radial.clip), stat,'uniformoutput',false)));
%	Sct = cvec(cell2mat(arrayfun(@(s) rvec(s.Sc.angular.clip), stat,'uniformoutput',false)));
	Sct = cvec(cell2mat(arrayfun(@(s) rvec(s.Sc.angular_resampled.pdf.clip), stat,'uniformoutput',false)));
	exclude      = cvec(arrayfun(@(x) x.exclude,stat));
	lcx          = cvec(arrayfun(@(x) 1./x.fc.x.clip,stat));
	lcr          = cvec(arrayfun(@(x) 1./x.fc.radial.clip,stat));
	isisotropic  = cvec(arrayfun(@(x) x.isisotropic,stat));
	ismodel      = cvec(arrayfun(@(x) x.ismodel,stat));
	is1d         = cvec(arrayfun(@(x) x.is1d,stat));
	isstochastic = cvec(arrayfun(@(x) x.isstochastic,stat));
	p_periodic   = cvec(arrayfun(@(x) x.p_periodic,stat));
	fhp          = cvec(arrayfun(@(x) x.fhp,stat));

	% relative domain size
	Le_rel_x     = Le_x./lcx;
	Le_rel_r     = Le_r./lcr;
	Le_rel       = Le_rel_x;
	Le_rel(isisotropic) = Le_rel_r(isisotropic);

	% characteristic wavelength
	lc              = lcx;
	lc(isisotropic) = lcr(isisotropic);

	% density maxima
	Sc1 = Scx;	
	Sc2 = Scy;
	fdx = isisotropic==1;
	Sc1(fdx) = Scr(fdx);
	Sc2(fdx) = Sct(fdx).*lc(fdx);

	% regularity
	regularityx   = Scx./lcx;
	regularityr   = Scr./lcr;
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
	file_C{fdx}
	%disp(fdx);
	
	% test for differences between natural and model generated patterns
	fdx0 =  ismodel & (~is1d) & (~isstochastic) & ~exclude;
	fdx1 = ~ismodel & (~is1d) & (~isstochastic) & ~exclude;
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
	tab = table;
	fdx = (exclude == 0);
	fprintf('N total %d\n',sum(fdx));
	id = 1;
	tab{id,1} = "Total";
	tab{id,2} = sum(fdx);
	%tab{id,3} = mean(p_periodic(fdx)<p_test);
	%tab{id,4} = median(Sc1(fdx)./lc(fdx));
	%tab{id,5} = median(Sc2(fdx)./lc(fdx));

	fdx = (~exclude) & (nL*lc < lc_) & (ismodel == 0);
	printf('Nature total: %d %0.2f %0.2f %0.2f\n',sum(fdx),mean(p_periodic(fdx)<p_test), median(Sc1(fdx)./lc(fdx)),median(Sc2(fdx)./lc(fdx)));
	id = 2;
	tab(id,1) = {"Nature total"};
	tab{id,2} = sum(fdx);
	tab{id,3} = round(mean(p_periodic(fdx)<p_test),2);
	tab{id,4} = round(median(Sc1(fdx)./lc(fdx)),2);
	tab{id,5} = round(median(Sc2(fdx)./lc(fdx)),2);
	fdx = ~exclude & (nL*lc < lc_) & (ismodel == 0) & (isisotropic == 0);
	printf('Nature aniso: %d %0.2f %0.2f %0.2f\n',sum(fdx),mean(p_periodic(fdx)<p_test), median(Scx(fdx)./lc(fdx)),median(Scy(fdx)./lc(fdx)));
	id = 3;
	tab{id,1} = {'Nature aniso'};
	tab{id,2} = sum(fdx);
	tab{id,3} = round(mean(p_periodic(fdx)<p_test),2);
	tab{id,4} = round(median(Sc1(fdx)./lc(fdx)),2);
	tab{id,5} = round(median(Sc2(fdx)./lc(fdx)),2);

	fdx = ~exclude & (nL*lc < lc_) & (ismodel == 0) & (isisotropic == 1);
	printf('Nature   iso: %d %0.2f %0.2f %0.2f\n',sum(fdx), mean(p_periodic(fdx)<p_test), median(Scr(fdx)./lc(fdx)),median(Sct(fdx)));
	id = 4;
	tab{id,1} = {'Nature iso'};
	tab{id,2} = sum(fdx);
	tab{id,3} = round(mean(p_periodic(fdx)<p_test),2);
	tab{id,4} = round(median(Sc1(fdx)./lc(fdx)),2);
	tab{id,5} = round(median(Sc2(fdx)./lc(fdx)),2);

	fdx = (exclude == 0) & ~exclude & (nL*lc < lc_) & (cvec(isstochastic) == 0) & (ismodel == 1) & (is1d == 0); 
	fprintf('2D Model Total: %d %0.2f %0.2f %0.2f\n',sum(fdx),mean(p_periodic(fdx)<p_test), median(Sc1(fdx)./lc(fdx)),median(Sc2(fdx)./lc(fdx)));
	id = 5;
	tab{id,1} = {'2D model deterministic total'};
	tab{id,2} = sum(fdx);
	tab{id,3} = round(mean(p_periodic(fdx)<p_test),2);
	tab{id,4} = round(median(Sc1(fdx)./lc(fdx)),2);
	tab{id,5} = round(median(Sc2(fdx)./lc(fdx)),2);


	fdx = (exclude == 0) & ~exclude & (nL*lc < lc_) & (cvec(isstochastic) == 0) & (ismodel == 1) & (isisotropic == 0) & (is1d == 0); 
	fprintf('2D Model Aniso: %d %0.2f %0.2f %0.2f\n',sum(fdx),mean(p_periodic(fdx)<p_test), median(Scx(fdx)./lc(fdx)),median(Scy(fdx)./lc(fdx)));
	id = 6;
	tab{id,1} = {'2D model aniso'};
	tab{id,2} = sum(fdx);
	tab{id,3} = round(mean(p_periodic(fdx)<p_test),2);
	tab{id,4} = round(median(Sc1(fdx)./lc(fdx)),2);
	tab{id,5} = round(median(Sc2(fdx)./lc(fdx)),2);

	fdx = (exclude == 0) & ~exclude & (nL*lc < lc_) & (cvec(isstochastic) == 0) & (ismodel == 1) & (isisotropic == 1);
	fprintf('2D model iso : %d %0.2f %0.2f %0.2f\n',sum(fdx), mean(p_periodic(fdx)<p_test), median(Scr(fdx)./lc(fdx)),median(Sct(fdx)));
	id = 7;
	tab{id,1} = {'2D model iso'};
	tab{id,2} = sum(fdx);
	tab{id,3} = round(mean(p_periodic(fdx)<p_test),2);
	tab{id,4} = round(median(Sc1(fdx)./lc(fdx)),2);
	tab{id,5} = round(median(Sc2(fdx)./lc(fdx)),2);

	fdx = (exclude == 0) & ~exclude & (nL*lc < lc_) & (cvec(isstochastic) == 1) & (ismodel == 1) & (isisotropic == 1);
	fprintf('2D Model stoch iso: %d %0.2f %0.2f %0.2f\n',sum(fdx), mean(p_periodic(fdx)<p_test),median(Scr(fdx)./lc(fdx)),median(Sct(fdx)));
	id = 8;
	tab{id,1} = {'2D model stoch iso'};
	tab{id,2} = sum(fdx);
	tab{id,3} = round(mean(p_periodic(fdx)<p_test),2);
	tab{id,4} = round(median(Sc1(fdx)./lc(fdx)),2);
	tab{id,5} = round(median(Sc2(fdx)./lc(fdx)),2);

	fdx = (exclude == 0) & ~exclude & (nL*lc < lc_) & (cvec(isstochastic) == 0) & (ismodel == 1) & (is1d == 1);
	fprintf('1D Model: %d %0.2f %0.2f\n',sum(fdx),mean(p_periodic(fdx)<p_test),median(Scx(fdx)./lc(fdx)));
	id = 9;
	tab{id,1} = {'1D model'};
	tab{id,2} = sum(fdx);
	tab{id,3} = round(mean(p_periodic(fdx)<p_test),2);
	tab{id,4} = round(median(Sc1(fdx)./lc(fdx)),2);
	tab{id,5} = round(median(Sc2(fdx)./lc(fdx)),2);
	tab.Properties.VariableNames = {'Group','N','$p<$0.05','$S_{c1}/\lambda_c$','$S_{c2}/\lambda_c$'};
	disp(tab);

	fid2 = fopen('mat/metastudy-table-2.tex','w');
	s = table2tex(tab);
	%fprintf(fid2,'\\begin{table}[H]
	%fprintf(fid2,'\\centering');
	fprintf(fid2,s);
	%fprintf(fid2,'\n\\end{table}\n');
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
		f = ['',dirname(f),'/crop/',basename(f)];
		subscript1 = ['xr'];
		subscript2 = {'y','\theta'};
		if (~exclude(idx))
		j = j+1;
		if (~is1d(idx))
			img2d = sprintf('	  \\includegraphics[height=0.23\\textwidth]{%s-density-2d-4-crop.pdf}',f); ...
			imgy  = sprintf('	  \\includegraphics[height=0.23\\textwidth]{%s-density-Sy-4-crop.pdf}',f); ...
			if (isisotropic(idx))
				regy_label = 'S_\theta';
				regy = stat(idx).Sc.radial.clip;
			else
				regy_label = 'S_y/\lambda_c';
				regy = stat(idx).Sc.y.clip/lc(idx);
			end
		else
			img2d = 'X (1d)';
			imgy  = 'X (1d)';
			regy_label = '-'
			regy      = NaN;
		end
		fprintf(fid,[...
		'\\begin{tabular}{cccc|c|c|c}\n' ...
		'\\multicolumn{6}{l}{%s %s %s}\n' ...
		'\\\\\\hline\n' ...
		'Density $S_{2d}$ & Density $S_%s$ & Density $S_%s$ & $S_{c%s}/\\lambda_c$ & $%s$ & $L_{eff}/\\lambda_c$ & $p$ \n' ...
		],...	
		num2str(idx),...
		figid,basename(d),...
		subscript1(isisotropic(idx)+1),...
		subscript2{isisotropic(idx)+1},...
		subscript1(isisotropic(idx)+1),...
		regy_label);
		fprintf(fid,[...	
		   '\\\\\\hline\n' ...
		'  %s\n' ...
		'& \\includegraphics[height=0.23\\textwidth]{%s-density-Sx-4-crop.pdf}\n' ...
		'& %s\n'], ...
		img2d,f,imgy);
	
		fprintf(fid,[ ...
		'& \\raisebox{11ex}{%0.2f}\n' ...
		'& \\raisebox{11ex}{%0.2f}\n' ...
		'& \\raisebox{11ex}{%0.1f}\n' ...
		'& \\raisebox{11ex}{%0.3f}\n' ...
		... '\\end{tblr}' 
		'\\end{tabular}\n\\\\' ...
		],...
		regularity(idx),regy,Le_rel(idx),p_periodic(idx));
	
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
		tab.regularity_y(j) = regy;
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

	end % if ~exclude
	
	end % for idx
	fclose(fid);
	writetable(tab,'mat/patterns-literature-stat.csv');

	% plot regularity	
	for xy=0:1
	for isiso_=0:1
	splitfigure([2,2],[30,1+isiso_+2*xy],fflag);
	cla();
	
	for ismodel_=0:1
	if (ismodel_)
		fdx =  ismodel & (~is1d) & (isisotropic == isiso_) & (~isstochastic) & ~exclude;
	else
		fdx = ~ismodel & (isisotropic == isiso_) & ~exclude;
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
			fdx_ = f>=0 & f<pi/2;
		else
			S = cvec(avgfun(Sr(fdx,:)));
			f = cvec(fi.x); %radial;
			fdx_ = f>=0;
		end
		S  = S/sum(mid(S(fdx_)).*diff(f(fdx_)));
		% nb: the circular should append the first to the last to f and Sa
		plot(f,S,'linewidth',1);
		if (xy)
			xlim([0,1]/4*2*pi);
			xlabel('Angle $\theta$','interpreter','latex');
			ylabel('Density $S_\theta$','interpreter','latex');
			set(gca,'xtick',[0,1/4,1/2]*pi,'xticklabel',{'0','\pi/4','\pi/2'});
		else
			xlim([0,2.5]);
			xlabel('Wavenumber $k_r / k_c$','interpreter','latex');
			ylabel('Density $S_r/\lambda_c$','interpreter','latex');
		end
	else
		if (xy)
			Sy(Sy<0) = 0;
			S = avgfun(Sy(fdx,:));
			f = fi.y;
		else
			% TODO what is going wrong here?
			Sx(Sx<0) = 0;
			S = avgfun(Sx(fdx,:));
			f = fi.x;
		end
		fdx = f>=0;
		df = (f(2)-f(1));
		S  = S/(sum(mid(S(fdx)))*df);
		plot(f,S,'linewidth',1);
		xlim([0,2.5]);
		if (xy)
		xlabel('Wavenumber $k_y/k_c$','interpreter','latex');
		ylabel('Density $S_y / \lambda_c$','interpreter','latex');
		else
			xlim([0,2.5]);
			xlabel('Wavenumber $k_x / k_c$','interpreter','latex');
			ylabel('Density $S_x/\lambda_c$','interpreter','latex');
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
	fdx = ismodel==0 & (~is1d) & (isisotropic == 0) & (~isstochastic) & ~exclude;
	c(1) = kendall_to_pearson(corr(cvec(Scx(fdx)),cvec(Scy(fdx))));

	figure(1003);
	clf();
	plot(log10(cvec(Scx(fdx))./lc(fdx)),log10(cvec(Scy(fdx)./lc(fdx))),'.')
	hold on
	fdx = ismodel & (~is1d) & (isisotropic == 0) & (~isstochastic) & ~exclude;
	c(2) = kendall_to_pearson(corr(cvec(Scx(fdx)),cvec(Scy(fdx))));
	plot(log10(cvec(Scx(fdx))./lc(fdx)),log10(cvec(Scy(fdx)./lc(fdx))),'.')
%	try
		fdx = ismodel==0 & (~is1d) & (isisotropic == 1) & (~isstochastic) & ~exclude;
		c(3) = kendall_to_pearson(corr(cvec(Scr(fdx))./lc(fdx),cvec(Sct(fdx))));
%	catch e
%	end

	figure(1004)
	clf();
	plot(log10(cvec(Scr(fdx))./lc(fdx)),cvec(Sct(fdx)),'.')
	hold on
	fdx = ismodel & (~is1d) & (isisotropic == 1) & (~isstochastic) & ~exclude;
	c(4) = kendall_to_pearson(corr(cvec(Scr(fdx))./lc(fdx),cvec(Sct(fdx))));
	plot(log10(cvec(Scr(fdx))./lc(fdx)),cvec(Sct(fdx)),'.')

	disp('Correlation:');
	disp(c)

	% plot correlation
	splitfigure([2,2],[100,1],fflag)
	cla();
	for idx=1:2
		plot(idx,c(idx),'*','color',cm(idx,:))
		hold on
	end
	ax = gca;
	set(ax(1),'ylim',[-0.2,1.09])
	ylabel(ax(1),'corr($S_x$,$S_y$)','interpreter','latex');
	xlim([0.5,2.5]);
	set(gca,'xtick',1:4,'xticklabel',{'nature','model','nature','model'},'xticklabelrot',45);
	grid on
	daspect([2,1,1])
%	yyaxis right


	% plot correlation
	splitfigure([2,2],[100,2],fflag)
	for idx=3:4
		plot(idx-2,c(idx),'*','color',cm(idx-2,:))
		hold on
	end
	ax=gca
	%set(ax,'ycolor','k');
	set(ax(1),'ylim',[-0.2,1.09])
	ylabel(ax,'corr($S_r$,$S_s$)','interpreter','latex');
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
		R_ = mean(Rx(fdx,:));
		plot(xi,R_,'linewidth',1);
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
		np(idx,1) = sum(fdx);
		q(:,idx) = quantile(regularity(fdx),[0.25,0.5,0.75]);
		qL_rel(:,idx)   = quantile(Le_rel(fdx),[0.25,0.5,0.75]);
	end
	da=[1.5   15.875    1];
	da=[1.5   10    1];
	da_=[1.5   0.375    1];
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
	ylabel('Regularity $S_{cx}/\lambda_c$','interpreter','latex');
	set(gca,'xtick',idp,'xticklabel',{'nature','2D-model','1D-model'},'xticklabelrot',45);
	set(gca,'ytick',2.^(-2:6));
	ylim([0.4,20]) 
	vline(2.5,'linestyle','--','color','k');
	ty = 1.4*32-20;
	%text(1-0.67*0-0.25,ty,'2D')
	%text(3-0.0,ty,'1D')
	daspect(da)
	
	% isotropic
	splitfigure([2,2],[10,2],fflag);
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
	ylabel('Regularity $S_{cr}/\lambda_c$','interpreter','latex');
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
		q(:,idx)  = quantile(cvec(Scy(fdx))./lc(fdx),[0.25,0.5,0.75]);
		fdx       = isisotropic & ~exclude & (ismodel == mo(idx)) & (isstochastic==sto(idx)) & (is1d == 0);
		qt(:,idx) = quantile(cvec(Sct(fdx)),[0.25,0.5,0.75]);
	end
	for idx=1:2
		h=errorbar(idx,q(2,idx),q(2,idx)-q(1,idx),q(3,idx)-q(2,idx),'*','color',cm(idx,:));
		hold on
	end
	%h=errorbar(3,q(2,3),q(2,3)-q(1,3),q(3,3)-q(2,3),'k*');
	xlim([0.5,2.5]);
	set(gca,'yscale','log');
	ylabel('Regularity $S_{cy}/\lambda_c$','interpreter','latex');
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
	ylabel('Regularity $S_{c\theta}$','interpreter','latex');
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
	
	splitfigure([2,2],[50,1],fflag);
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
		R(~isisotropic,:) = Rx(~isisotropic,:);
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
	set(gca,'colororder',meta.colormap);
	
	if (0)
		S1d_(~isfinite(S1d_))=0;
		m=0;
		R1d = real(ifft([S1d_; zeros(m,size(S1d_,2)); flipud(S1d_(2:end,:))]));
		R1d = R1d./R1d(1,:);
		df=f1d(2)-f1d(1);
		L = 1./df;
		x=linspace(0,L,length(R1d))';
	end
	
	% plot autocorrelation
	splitfigure([2,2],[50,2],fflag);
	cla
if (0)
	plot(xi,R1d_./R1d_(1,:),'linewidth',1);
end
	xlabel('Lag Distance x/\lambda_c');
	ylabel('Autocorrelation R');
	xlim([0,2.25])
	set(gca,'colororder',meta.colormap);
	
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
		R(~isisotropic,:) = Rx(~isisotropic,:);
		R1d_(:,idx) = nanmean(R(fdx,:));
end	
	end
	S1d_ = S1d_./(sum(S1d_)*(f1d(2)-f1d(1)));
	
	plot(f1d,S1d_,'linewidth',1);
	xlabel('Wavenumber k/k_c');
	ylabel('Density S_c/\lambda_c');
	legend('real','model','model 1d');
	xlim([0,2.5])
	set(gca,'colororder',meta.colormap);
	
	fdx = ~exclude & ~exclude;
	f = file_C(fdx);
	f = cellfun(@dirname,f,'uniformoutput',false)';
	f =unique(f)


	if (meta.pflag)
		ps = 3.5;
		ps_ =4;
		pdfprint(101,'img/metastudy-regularity-Scx.pdf',ps_);
		pdfprint(102,'img/metastudy-regularity-Scr.pdf',ps_);
		pdfprint(103,'img/metastudy-regularity-Scy.pdf',ps_);
		pdfprint(104,'img/metastudy-regularity-Sct.pdf',ps_);
		pdfprint(1001,'img/metastudy-regularity-correlation-Sx-Sy.pdf',ps_);
		pdfprint(1002,'img/metastudy-regularity-correlation-Sr-St.pdf',ps_);
	if (0)
		figure(103); axis square
		pdfprint(103,'img/density-literature.pdf',ps);
		pdfprint(104,'img/autocorrelation-literature.pdf',ps);
	
		pdfprint(201,'img/regularity-literature-simple.pdf',ps);
		pdfprint(203,'img/density-literature-simple.pdf',ps);
	%	pdfprint(204,'img/autocorrelation-literature-simple.pdf',ps);
	end
		pdfprint(301,'img/metastudy-density-Sx.pdf',ps);
		pdfprint(302,'img/metastudy-density-S-radial.pdf',ps);
		pdfprint(303,'img/metastudy-density-Sy.pdf',ps);
		pdfprint(304,'img/metastudy-density-S-angular.pdf',ps);
	if (0)
		pdfprint(401,'img/metastudy-autocorrelation-Rx.pdf',ps);
		pdfprint(402,'img/metastudy-autocorrelation-Rr.pdf',ps);
	end
		
	end % if pflag
	
end % pattern_metastudy_plot


