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
% batch script for reproducing the analysis, figures and tables
%
% for copyright reasons, the repository only includes images appearing in the
% accompanying publication, the images extracted from the references are not
% included

%	fetch_dependcies();

	addpath(['./lib/auxiliar']);
	addpath_recursive('./lib');

	mkdir('mat/');
	mkdir('img/');

	meta = pattern_metadata();

	% set to true to save fitures to files 
	pflag      = false;
	meta.pflag = pflag;

	pattern_analysis_minimum_working_example();
	
	pattern_observed_plot_2d(meta);

	pattern_schematic_plot(meta);

	pattern_regularity_sweep();

	pattern_metastudy_analyze([],[],meta.pflag);
	pattern_metastudy_plot(meta);

	experiment_density_averaging();

	experiment_regularity_vs_p_value();

	experiment_regularity_estimate_bias();

