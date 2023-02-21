% Wed  8 Feb 09:36:33 CET 2023
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
% minimum working example for analyzing the regularity of spatial patterns
%

	img_C = { 'patterns/2a_pattern_striped_11.33051_28.35545_model_0.png', ...
		  'patterns/2b_pattern_2d_+11.53386_+027.92788_model_0.png' ...
		};

	for idx=1:2

	% create spatial pattern analysis object
	sp = Spatial_Pattern();                                         
	
	% load image
        sp.imread(img_C{idx}); 

	% analyze
	sp.analyze_grid();

	% tabulate analysis results
	sp.report();

	figure(idx);
	clf();

	% plot pattern
	subplot(2,3,1);
	sp.plot('b');
	title('Pattern in real space');

	% plot 2D periodogram, low-frequency components suppressed
	subplot(2,3,2);
	sp.plot('S.clip');
	title('2D-Periodogram');

	% plot 2D density
	subplot(2,3,3);
	sp.plot('S.bar');
	title('2D-Density');

	if (sp.stat.isisotropic)
		% plot 1D density along primary axis
		subplot(2,3,4);
		sp.plot('S.radial.hat');
		hold on
		sp.plot('S.radial.clip');
		legend('complete','low-frequencies suppressed')
		title('Density along major axis');
		ylim([0,1])
	
		% plot 1D density along secondary axis
		subplot(2,3,5);
		sp.plot('S.rot.angular.bar');
		title('Density along secondary axis');
		ylim([0,1])
	else
	
		% plot 1D density along primary axis
		subplot(2,3,4);
		sp.plot('S.rot.x.hat');
		hold on
		sp.plot('S.rot.x.clip');
		legend('complete','low-frequencies suppressed')
		title('Density along major axis');
	
		% plot 1D density along secondary axis
		subplot(2,3,5);
		sp.plot('S.rot.y.clip');
		title('Density along secondary axis');
	end % else of if isstochastic

	end % for idx

