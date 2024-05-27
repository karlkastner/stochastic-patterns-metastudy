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
%% analyze the patterns in the metastudy
%%
%% attributes encoded in file name:
%% <folder>/<number>_
%%
%% folder : title,first author,year of respective publication
%% number : number of figure and subplot in the publication, e.g. 1a
%% model  : {0 = nature or physical experiment,
%%           1 = generated by numerical model}
%% 1d     : {0 = 1d-pattern (only for models),
%%           1 : 2d-pattern}
%% stoch  : {0 : model without stochastic components,
%%           1 : model with stochastic components}
%% exclude : {0 : include file in the analysis,
%%            1 : exclude file in the analysis}
%%
function [stat,file_C,Si,Ri] = pattern_metastudy_analyze(id,meta)
	if (nargin()<2)
		meta = pattern_metastudy_metadata();
	end

	pflag = meta.pflag;
	fflag = pflag;

	% number of subplots
	nsubplot = 5; 

	% number of values in the resamples densities (for plotting)	
	ni.x     = 400;
	ni.y     = 401;

	% resampled coordinate axes
	xi         = (0:(1/50):5)';
	fi.radial  = linspace(0,5,ni.x)';
	fi.x       = linspace(0,5,ni.x)';
	fi.y       = linspace(0,3.5,ni.y)';
	fi.angular = innerspace(-pi/2,pi/2,ni.y)';

	% patterns from the accompanying publication
	cmd = 'ls patterns/[0-9]*.*g -1'; % png and jpg patterns/[0-9]*.jpg -1';
	[~,ret_str] = system(cmd);
	file_C_ = strsplit(ret_str,'\n');

	% patterns from the literature
	cmd = 'ls patterns/metastudy/*/[0-9]*.png -1';
	[~,ret_str] = system(cmd);
	file_C = strsplit(ret_str,'\n');

	% concatenate
	file_C = [file_C_,file_C];

	% remove the empty strings due to last carriage return
	file_C = file_C(~cellfun(@isempty,file_C));

	% default : select all files for processing	
	lf = length(file_C);
	if (nargin()<1||isempty(id))
		id=(1:length(file_C));
	else
		id(id > lf) = [];
	end
	
	% allocate output variables
	clear stat Si Ri;
	
	for idx=id
		printf('Pattern: %d\n',idx);
	
		% figure and subplot index
		plotid = mod(idx-1,nsubplot)+1;
		if (~isempty(file_C{idx}))
	
		% parse properties from file name
		% if file was model generated
		model_str=regexprep(file_C{idx},'.*model_([0-9]*).*','$1');
		ismodel = str2double(model_str);

		% if file is to be excluded, usually bc lenght-scale is not clearly separated from 0
		if (regexp(file_C{idx},'exclude_[0-9]'))
			exclude_str=regexprep(file_C{idx},'.*exclude_([0-9]*).*','$1');
			exclude = str2double(exclude_str);
		else
			exclude = 0;
		end
		% if the pattern was generated by a stochastic model
		if (regexp(file_C{idx},'stoch_[0-9]'))
			isstoch_str=regexprep(file_C{idx},'.*stoch_([0-9]*).*','$1');
			isstochastic = str2double(isstoch_str);
		else
			isstochastic = 0;
		end
		if (regexp(file_C{idx},'sdf_[0-9]'))
			issdf_str=regexprep(file_C{idx},'.*sdf_([0-9]*).*','$1');
			issdf = str2double(issdf_str);
		else
			issdf = 1;
		end
		% if the pattern was generated by a 1D-model
		if (regexp(file_C{idx},'1d_[0-9]'))
			%~isempty(d1_str))
			is1d_str = regexprep(file_C{idx},'.*1d_([0-9]*).*','$1');
			is1d = str2double(is1d_str);
		else
			is1d = 0;
		end

		% reset random number generator for reproducibility of test
		rng(0);
		sp = Spatial_Pattern();
		% analyze pattern
		sp.imread(file_C{idx});	
		sp.analyze_grid();
		% resample densities to common axes
		[Si(idx,1),Ri(idx,1)] = sp.resample_functions(xi,fi);
		% fit parameteric densities
		try
			sp.fit_parametric_densities();
		catch e
			disp(e);
			sp.stat.fit = struct();
		end

		% pattern properties
		if (sp.stat.isisotropic)
			%Sc = sp.stat.Sc.radial.hp;
			lc = 1./sp.stat.fc.radial.hp;	
			Le_rel = sp.stat.L_eff.r ./ lc;
		else
			%Sc = sp.stat.Sc.x.hp;
			lc = 1./sp.stat.fc.x.hp;	
			Le_rel = sp.stat.L_eff.x ./ lc;
		end % else of isiso

	if (meta.dflag)
		% plot annotation	
		title_C = { sprintf('%d p_{iso}=%1.2f, S_c/l_c=%1.2f p=%0.2f nf=%d',idx,...
									   sp.stat.p_isotropic,...
									   sp.stat.regularity,...
									   sp.stat.p_periodic,...
									   sp.stat.nf),...
			    sprintf('Le/lc = %0.1f mod=%d excl=%d',Le_rel,ismodel,exclude) };
		figid = floor((idx-1)/nsubplot)+1;
	
		% plot masked pattern
		splitfigure([5,nsubplot],[figid,plotid],fflag,[],100,false,[],'Visible',meta.visible);
		cla();
		sp.plot('b');
		title(title_C{:});
	
		% plot periodogram
		splitfigure([5,nsubplot],[figid,plotid+nsubplot],fflag,[],100,false,[],'Visible',meta.visible);
		cla();
		if (~isfinite(lc))
			lc_ = 1./sp.stat.q.fr.p50;
		else
			lc_ = lc;
		end
		sp.plot('S.rot.hat');
		hold on
		% mark the region tested for periodicity
		contour(fftshift(sp.f.x*lc_),fftshift(sp.f.y*lc_),fftshift(sp.msk.rot.f)',[-1,0.5],'r')
		contour(fftshift(sp.f.x*lc_),fftshift(sp.f.y*lc_),fftshift(sp.f.rr<sp.stat.fhp)',[-1,0.5],'g','linewidth',1)
		
				
		splitfigure([5,nsubplot],[figid,plotid+2*nsubplot],fflag,[],100,false,[],'Visible',meta.visible);
		sp.plot('S.rot.bar');
		if (~fflag)
		hold on
		contour(fftshift(sp.f.x*lc_),fftshift(sp.f.y*lc_),fftshift(sp.msk.rot.f)',[-1,0.5],'r');
		end
		if (~fflag)
			colormap(parula);
		end
		
		% plot one-dimensional density
		splitfigure([5,nsubplot],[figid,plotid+3*nsubplot],fflag,[],100,false,[],'Visible',meta.visible);
		cla();
		if (sp.stat.isisotropic)
			sp.plot('S.radial.hp');
			hold on
			sp.plot('S.radial.hat');
			if (sp.stat.fc.radial.hp>0)
			ylim([0,1.05*sp.stat.Sc.radial.hp*sp.stat.fc.radial.hp]);
			end
			vline(sp.stat.fhp/sp.stat.fc.radial.hp);
		else
			sp.plot('S.rot.x.hp');
			hold on
			sp.plot('S.rot.x.hat');
			if (sp.stat.fc.x.hp>0)
			ylim([0,1.05*sp.stat.Sc.x.hp*sp.stat.fc.x.hp]);
			end
			vline(sp.stat.fhp./sp.stat.fc.x.hp);
		end % else of isisotropic
		vline(sp.stat.fhp*lc,'color','k','linestyle','--');
		if (fflag)
			axis square
		end
	
		% plot one-dimensional density on secondary axis,
		% i.e. parallel to bands or angular	
		splitfigure([5,nsubplot],[figid,plotid+4*nsubplot],fflag,[],100,false,[],'Visible',meta.visible);
		cla()
		if (~is1d)
		if (sp.stat.isisotropic)
			sp.plot('S.rot.angular.hp');
			hold on
			yl = ylim;
			ylim([0,max(1,yl(2))]);
		else
			sp.plot('S.rot.y.hp');
		end
		end % else of is1d
		if (fflag)
			axis square
		end

		if (pflag)
			ps = 4;
			base = basename(file_C{idx});
			dir  = dirname(file_C{idx});
			dir_C = strsplit(dir,filesep);
			if (length(dir_C)>1)
				dir  = [dir_C{end-1},'/',dir_C{end}];
			else
				dir = dir_C{1};
			end
			imgbase = ['img/metastudy-result/',dir,'/',base(1:end-4)]
			mkdir(imgbase);
			pdfprint(100*figid+plotid+1*nsubplot,[imgbase,'-periodogram-2d.pdf'],ps);
			pdfprint(100*figid+plotid+2*nsubplot,[imgbase,'-density-2d.pdf'],ps);
			pdfprint(100*figid+plotid+3*nsubplot,[imgbase,'-density-Sx.pdf'],ps);
			pdfprint(100*figid+plotid+4*nsubplot,[imgbase,'-density-Sy.pdf'],ps);
			close all;
		end % if pflag
	end  % if dflag
	
	end % if not isempty file_C{idx}

	% store analysis results
	sp.stat.exclude         = exclude;
	sp.stat.is1d            = is1d;
	sp.stat.ismodel         = ismodel;
	sp.stat.isstochastic    = isstochastic;
	sp.stat.issdf           = issdf;
	stat(idx,1)             = sp.stat;

	end % for idx
	
	if (meta.saveflag)
		save(meta.filename.metastudy,'file_C','stat','xi','fi','Si','Ri');
	end
end % function pattern_metastudy_analyze

