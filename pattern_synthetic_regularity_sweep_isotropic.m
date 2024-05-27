% Thu 23 Feb 15:25:26 CET 2023
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
%% synthesize a spatial isotropic pattern with gradually varying regularity
%
if (1)
if (~exist('pflag','var'))
	pflag = 0;
end
ps = 1.5;

% characteristic wavelength of patern
lc = 1;
% length of domain, determining spectral resolution
L = 40*lc;
m = 40;
% length of grid cells, determining spatial resolution
% the spectral and spatial resolution with respect to the pattern wavelength
% is kept identical here, but this is not necessary
dx = lc/m;
% number of grid points in one dimension
n = L/dx;
% phase noise scale, determining the regularity of the synthesized patterns
s   = flipud(2.^(-2:0.5:3.5)');
% seed of random number generator
seed = 1;

rng_ = 1;
% output file name
filename=sprintf('mat/synthetic-isotropic-pattern-rng-%d-L-%d-m-%d.mat',rng_,L,m);
filename
if (exist(filename,'file'))
	disp('Loading file');
	load(filename);
else
	bb  = [];
	Scr = [];
	Sca = [];
	
	% generate several patterns with the same background noise but different regularity
	for idx=1:length(s)
		disp(idx)
		% reset random number generated so that all patterns are generated with the same noise
		rng(rng_);
		% generate hexagonal pattern with phase noise
		%[b,x,y]=generate_isotropic_pattern(1/lc,n,L,0,0,s(idx),0.5,1);
		alpha = 0;
		[b,x,y]=generate_isotropic_pattern(1/lc,n,L,alpha,[],[],0,s(idx));
		bb(:,:,idx) = b;
		% 2d-periodogram
		S = abs(fft2(b-mean(b(:)))).^2;
		% spectral density
		Sbar=trifilt2(S,3);
		% radial density
		[Sr,fr] = periodogram_radial(S,[L,L]);
		% determine regularity
		fdx=find(fr==1);
		Scr(idx,1) = Sr.normalized(fdx);
		Scr(idx,2) = max(Sr.normalized);
		% angular density
		[Sa,angle] = periodogram_angular(S,[L,L]);
		Sca(idx,1) = interp1(angle,Sa,0,'linear');
		Sca(idx,2) = max(Sa);
		
		if (0)
			figure(1);
			subplot(4,6,idx);
			imagesc(b>0.7);
			axis equal;
			axis tight;
			colormap(flipud(gray));
			figure(2);
			subplot(4,6,idx);
			imagesc(fftshift(Sbar));
			axis equal;
			axis tight;
			%plot(fr,Sr.normalized);
			%vline(fr(fdx))
			title([max(Sr.normalized),Sr.normalized(fdx)]);
			figure(3);
			subplot(4,6,idx);
			plot(fr,Sr.normalized);
			figure(4);
			subplot(4,6,idx)
			plot(angle,Sa);
		end
	end % for idx
	
	% bb contains several patterns with similar arrangement of patches,
	% but different regularity, these are now combined into one pattern
	% where the regularity gradually increases from left to right
	disp('Combining patterns with different regularity into one pattern with gradually varying regularity');
	Sc_ = sort(Scr(:,1));
	iSc = logspace(log10(0.125),log10(16),n);
	bi = [];
	for idx=1:n
		for jdx=1:n
			bi(idx,jdx) = interp1(Sc_,squeeze(bb(idx,jdx,:)),iSc(jdx),'linear');
		end
	end
	% save the pattern
	save(filename,'x','y','Scr','Sca','s','bb','bi');
end % else of if exist filename
	iSc = logspace(log10(0.125),log10(16),n);

p_periodic = [];
for idx=1:size(bb,3)
	nf_test = 3;
	bmsk = [];
	fmsk =[];
	[isperiodic(idx), p_periodic(idx), stati, out] = periodogram_test_periodicity_2d(...
								bb(:,:,idx), [L,L], nf_test, bmsk, fmsk); 
end
p_periodic
Sc_ = sort(Scr(:,1));                                                   
p_periodic_i = interp1(Sc_,p_periodic,iSc,'linear');
pdx=find(p_periodic_i<=0.05,1,'first');


figure(6)
clf;
plot(1./s,[Scr,Sca]/lc);
ylabel('Phase noise scale s');
ylabel('Regularity $Sc/lc$','interpreter','latex');

% plot threshholded pattern
figure(5);
imagesc(x,y,bi>0.6);
axis equal;
axis tight;
Sct = 2.^(-3:4);
xt  = interp1(log2(iSc),x,log2(Sct),'linear');
set(gca,'xtick',xt,'xticklabel',num2str(cvec(Sct)))
colormap(flipud(gray));
xlabel('Regularity $S_{cr}/\lambda_c$','interpreter','latex');
ylabel('Distance $y/\lambda_c$','interpreter','latex');
colormap(flipud(colormap_vegetation()));
%vline(iSc(pdx)/lc,'linewidth',1.5,'color','r','linestyle','-')

dat = load('mat/patterns-metastudy.mat');                     
%lc_ = cvec(1./arrayfun(@(x) x.fc.radial.hp,dat.stat));                   
ismodel = cvec(arrayfun(@(x) x.ismodel,dat.stat));                           
isisotropic = cvec(arrayfun(@(x) x.isisotropic,dat.stat));                      
exclude = cvec(arrayfun(@(x) x.exclude,dat.stat));                              
Scr_ =    cvec(arrayfun(@(x) x.Sc.radial.hp,dat.stat));                               
lcr_ = 1./cvec(arrayfun(@(x) x.fc.radial.hp,dat.stat));                               
%Scy_ = cvec(arrayfun(@(x) x.Sc.y.clip,dat.stat));  
col = 'br';
end
for idx=1:2
	fdx = ismodel == (idx-1) & (isisotropic == 1) & exclude == 0;
	Scr_lcr_q = quantile(cvec(Scr_(fdx)./lcr_(fdx)),[0.25,0.5,0.75])
	% interpolate to x
	xi = interp1(iSc/lc,x,Scr_lcr_q,'linear');

	vline(xi(1),'color',col(idx),'linewidth',1.5);
	vline(xi(3),'color',col(idx),'linewidth',1.5);
end

% save the figure, cropping the y-axis differently
if (pflag)
Lp = 10:10:L;
for idx=1:length(Lp)
	figure(5);
	ylim(Lp(idx)/2*[-1,1]);
	pdfprint(5,sprintf('img/synthetic-isotropic-pattern-L-%d-m-%d-Lp-%d.pdf',L,m,Lp(idx)),ps);
end
end


