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
%	dependencies_fetch(meta.url,meta.filename.dependencies);

	% add libraries to path
	addpath_recursive('./lib');

	% set to true to save figures to files 
	pflag      = false;
	meta.pflag = pflag;

	pattern_analysis_minimum_working_example();

	% insets for figures of a schematic anisotropic and isotropic pattern
	plot_schematic_patterns_insets();

	% Figure  1	
	pattern_observed_plot_2d(meta);

	% Figure 2
	plot_schematic_periodic_vs_stochastic();

	% Figure 3
	plot_schematic_density_decomposition_2d;

	% Figure 4
	experiment_spectral_density_scaling();

	% Figure 5 a-b and Figure 6 a-b
	plot_schematic_density_decomposition_1d();

	% Figure 5 c
	plot_anisotropic_density_2d();

	% Figure 5 d
	pattern_synthetic_regularity_sweep_anisotropic();

	% Figure 6 c
	pattern_synthetic_regularity_sweep_isotropic();

	% the metastudy requires images with patterns from the references
	pattern_metastudy_analyze([],meta);

	% Figure 7 and Figure 8
	pattern_metastudy_plot(meta);

	% Figure SI  1
	experiment_area_aspect_ratio();

	% Figure SI  2 and Figure SI 3	
	experiment_regularity_estimate_uncertainty_components();

	% Figure SI  4
	experiment_regularity_estimate_uncertainty_varying_extent();

	% Figure SI  5
	experiment_density_averaging();

	% Figure SI  6
	experiment_regularity_vs_p_value();

	% Figure SI  7
	example_patterns_several_distributions();

	% Figure SI  8
	experiment_regularity_measure_convergence();

	% Figure SI  9 ai-aiii
	experiment_regularity_measure_finite_spatial_extent();	

	% Figure SI  9 bi-biii
	experiment_regularity_measure_finite_sampling_interval();

	% Figure SI 10
	experiment_regularity_estimate_uncertainty_measure_comparison();

