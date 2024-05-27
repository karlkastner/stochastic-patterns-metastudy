% Mon 31 May 20:20:46 CEST 2021
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
%% metadata
%
function meta = pattern_metastudy_metadata()
	meta.pflag = 0;

	meta.saveflag = true;
	meta.dflag = 1;
	meta.visible = 'off';


	meta.p_test = 0.05;
	
	meta.url                   = 'https://github.com/karlkastner/';

	meta.filename.dependencies = 'dependencies.csv';
	meta.filename.profile      = 'mat/profiling-information.mat';
	meta.filename.metastudy    = 'mat/patterns-metastudy.mat';

end

