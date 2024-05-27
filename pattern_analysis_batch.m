% 2021-11-26 19:40:22.060116517 +0100
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
%% batch script for reproducing the analysis, figures and tables
%
	meta = pattern_metastudy_metadata();
	
	% create library and output folder
	mkdir('./lib/');
	mkdir('./mat/');
	mkdir('./img/');
	mkdir('./lib/auxiliar/');
	addpath('./lib/auxiliar');

	% fetch the script for fetching library files
	%cmd = sprintf(['svn export %s/auxiliar/trunk/dependencies_fetch.m ./lib/auxiliar/'],meta.url);
	%system(cmd);
	url  = 'https://raw.githubusercontent.com/karlkastner/auxiliar/master/dependencies_fetch.m';
	dest = './lib/auxiliar/dependencies_fetch.m';
	urlwrite(url,dest);

	% fetch library files
	% dependencies_determine(meta.filename.dependencies,meta.filename.profile,{'pattern_analysis_batch','pdfprint'});
	dependencies_fetch(meta.url,meta.filename.dependencies);

	% add libraries to path
	addpath_recursive('./lib');

	% set to true to save figures to files 
	pflag      = false;
	meta.pflag = pflag;

	pattern_analysis_minimum_working_example();
	
	pattern_observed_plot_2d(meta);

	pattern_synthetic_plot(meta);

	pattern_anisotropic_regularity_sweep();

	pattern_isotropic_regularity_sweep();

	experiment_density_averaging();

	experiment_regularity_vs_p_value();

	experiment_regularity_estimate_bias();

	% the metastudy requires images with patterns from the references
	pattern_metastudy_analyze([],[],meta.pflag);
	pattern_metastudy_plot(meta);

