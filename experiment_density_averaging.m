% 2023-02-10 16:50:21.623554632 +0100
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
% demonstration, that averaging densities with the same distribution
% but different regularity results in a density that is more pointed
% and has heavier tales than the underlying distribution
%

% spatial extend
L = 50;
% number of samples in one dimension
n = 50^2;
% frequencies
fx = linspace(0,L,n)';
% number of densities to average
m = 1000;

% model the regularity log-normally distributed
% and systematically (not randomly) sample between 0 and 1 from the cdf
p = innerspace(0,1,m);
fc0 = 1;
Sc0 = 1;
s = 1;
[a,b] = logn_moment2param(Sc0,s*Sc0);
Sc = logninv(p,a,b);
cv = std(Sc)/mean(Sc)

% the length scale can also be varied, after normalization, this results
% however again in a variation of Sc
s = sqrt(eps);
[a_,b_] = logn_moment2param(fc0,s*fc0);
fc = logninv(p,a_,b_);
p = [];
S = [];
df = 1/L;
Sflat = 1/(n*df);

% generate the distributions with different regularity
p = 1-sqrt(eps);
for idx=1:length(Sc)
	[a,b] = logn_mode2param(fc(idx),Sc(idx));
	S(:,idx) = p*lognpdf(fx*fc(idx),a,b)*fc(idx)+(1-p)*Sflat;
%	[a,b] = gamma_mode2par(fc,Sc(idx));
%	S(:,idx) = gampdf(fx,a,b);
end
% normalize
S = S./(sum(S)*df);
% average
S_ = [mean(S,2)];
S_ = S_./(sum(S_)*df);
% distributions for comparison which are not averaged
[a,b] = logn_mode2param(fc0,Sc0*[0.5,1,2]);
S1 = lognpdf(fx,a,b);
ls = {'--','-','-.'};

% display
figure(1)
clf()
for idx=1:3
	plot(fx,S1(:,idx),['k',ls{idx}],'linewidth',1);
hold on
end
plot(fx,S_,'r','linewidth',1);
xlim([0,2.5]);
xlabel('Wavenumber k/k_c');
ylabel('Density S/\lambda_c');
lh=legend('1/2','1','2','avg')
title(lh,'Regularity S_c/\lambda_c')

if (exist('pflag','var') && pflag)
	pdfprint(1,'img/density-averaging.pdf',3.5);
end

