% Sat 28 Jan 23:05:35 CET 2023
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
%% illustrate the visual appearance of a pattern depending on the regularity
%% on a synthesized a pattern where the regularity varies
%

if (~exist('pflag','var'))
	pflag = 0;
end
mode     = 'relxy';
pdf_str  = 'logn';
pdfy_str = 'exp';
ps = 1.5;
% mode='rely';

% distribution along the primary
pdf = @lognpdf;
pdf_mode2par = @logn_mode2par;

% distribution along the secondary axis 
pdfy = @exppdf;
pdfy_mode2par = @exppdf_max2par;

% spatial extend
L  = [60,10];
% characteristic wavelenght
lc = 1;
% characteristic frequency
fc = 1/lc;
% spatial discretization
dx = lc/10;
% number of grid point
n  = L/dx;
% reset random number generator (for reproducibility)
rng_ = 0;
rng(rng_)

switch (mode)
case {'equal'}
	Scx = 2*logspace(-1,1,n(1))';
	Scy = logspace(-1,1,n(1));
case {'independent'}
	Scx = 2*logspace(-1,1,n(1))';
	L(2) = L(1);
	n(2) = n(1);
	Scy = logspace(-1,1,n(2));
case {'rely'}
	L(2) = 30;
	n(2) = m*L(2);
	Scy_rel = logspace(-1,1,n(2));
	Scy = cvec(Scx)*Scy_rel;
case {'relxy'}
	ScxScy_lim       = [0.125,16];
	Scy_div_Scx_lim  = [2.^-1.5,2.^1.5];
	L = [60,30];
	L(2) = L(1)*range(log10(Scy_div_Scx_lim))/range(log10(ScxScy_lim));
	m = 10;
	n  = round(m*L);
	ScxScy = logspace(log10(ScxScy_lim(1)),log10(ScxScy_lim(2)),n(1))';
	Scy_div_Scx = logspace(log10(Scy_div_Scx_lim(1)),log10(Scy_div_Scx_lim(2)),n(2))';
	Scx = ScxScy./rvec(Scy_div_Scx);
	Scy = ScxScy.*rvec(Scy_div_Scx);

end

% output file
%f_str = sprintf('pattern-sweep-%s-%s-%s-Scx-%f-%f-Scy-%f-$f-L-%f-%f-n-%d-%d',pdf_str,pdfy_str,mode,min(Scx(:)),max(Scx(:)),min(Scy(:)),max(Scy(:)),L,n);
filename=sprintf('mat/synthetic-anisotropic-pattern-rng-%d-Scx-%f-%f-Scy-%f-%f-L-%d-%d-dx-%d.mat',rng_,min(Scx(:)),max(Scx(:)),min(Scy(:)),max(Scy(:)),L(1),L(2),dx(1));
%f_str = ['mat/',f_str,'.mat'];

if (exist(filename,'file'))
	disp('Loading file');
	load(filename);
else

% axes in real space
x = innerspace(0,L(1),n(1));
y = innerspace(0,L(2),n(2));

% axes in frequency space
[fx,fy,fr] = fourier_axis_2d(L,n);

e  = randn(n);
fe = fft2(e);

% vary y-regularity on y-axis
%-> vary Scy on y-axis

% compute density parameters a and b from characteristic wavelength 1/fc and regularity Sc*fc 
b  = zeros(n);
br = [];
Sx_ = zeros(n(2),n(2));
Sy_ = zeros(n(1),n(1));
p = struct();
flag = false;
for idx=1:n(1)
	for jdx=1:n(2)
		disp([idx/n(1),jdx/n(2)])
	if (1)
		[p.x(idx,jdx,1),p.x(idx,jdx,2)] = pdf_mode2par(fc,Scx(idx,jdx)); 
		if (~strcmp(pdfy_str,'exp'))
			[p.y(idx,jdx,1),p.y(idx,jdx,2)] = pdfy_mode2par(0,Scy(idx,jdx)); 
		else
			[p.y(idx,jdx,1)] = pdfy_mode2par(Scy(idx,jdx)); 
		end
	else
		[p.x(idx,jdx,1),p.x(idx,jdx,2)] = pdf_mode2par(fc,Scx(idx),[p.x(idx-1,jdx,1),p.x(idx-1,jdx,2)],flag);
		[p.y(idx,jdx,1),p.y(idx,jdx,2)] = pdf_mode2par(0,Scy(idx),[p.y(idx-1,jdx,1),p.y(idx-1,jdx,2)],flag);
	end
	end % for jdx
end % for idx

switch (mode)
	case {'equal'}
	if (equal)
	%for idx=1:n(1)
		% spectral density
		Sx = gampdf_man(abs(fx),p(idx,1),p(idx,2));
		Sy = gampdf_man(abs(fy),q(idx,1),q(idx,2));
		% transfer function
		Tx = sqrt(Sx);
		Ty = sqrt(Sy);
		% 2D transfer function
		T  = (cvec(Ty)*rvec(Tx));
		bi = ifft2(T.*fe);
		b(:,idx) = real(bi(:,idx)); 
	
		Sr = gampdf_man(fr,ap,bp);
		Tr = sqrt(Sr);
		bi = ifft2(Tr.*fe);
		br(:,idx) = real(bi(:,idx));
	end
	case {'independent'}
	br = zeros(n);
	for idx=1:n(1)
	disp(idx/n(1));
	for jdx=1:n(2)
		% spectral density
		Sx = gampdf(abs(fx),p(idx,1),p(idx,2));
		Sy = gampdf(abs(fy),q(jdx,1),q(jdx,2));
		% there is a bug in matlab
		Sy(1) = 2*Sy(2)-Sy(3);
	
		% transfer function
		Tx = sqrt(Sx);
		Ty = sqrt(Sy);
		% 2D transfer function
		T  = (cvec(Ty)*rvec(Tx));
		T = fftshift(T);
		T = imrotate(T,-45,'crop','bilinear');
		T = ifftshift(T);
		bi = ifft2(T.*fe);
		b(idx,jdx) = real(bi(idx,jdx)); 
	end % for idx
	%if (0)
		Sr = gampdf_man(fr,p(idx,1),p(jdx,2));
		Tr = sqrt(Sr);
		bi = ifft2(Tr.*fe);
		br(idx,:) = real(bi(idx,:));
	%end
	end % for jdx
	case {'relxy','rely'}

	isperiodic = [];	
	p_periodic = [];
	for idx=1:n(1)
		disp(idx/n(1));
		for jdx=1:n(2)
			% spectral density
			Sx = pdf(abs(fx),p.x(idx,jdx,1),p.x(idx,jdx,2));
			if (~strcmp(pdfy_str,'exp'))
				Sy = pdfy(abs(fy),p.y(idx,jdx,1),p.y(idx,jdx,2));
			else
				Sy = pdfy(abs(fy),p.y(idx,jdx,1));
			end
			% there is a bug in matlab at gammapdf at 0
			% Sy(1) = 2*Sy(2)-Sy(3);
		
			% transfer function
			Tx = sqrt(Sx);
			Ty = sqrt(Sy);
			% 2D transfer function
			T  = (cvec(Tx)*rvec(Ty));
			bi = ifft2(T.*fe);
			nf_test = 3;
			bmsk = [];
			fmsk =[];
			[isperiodic(idx,jdx), p_periodic(idx,jdx), stati, out] = periodogram_test_periodicity_2d(...
								bi, L, nf_test, bmsk, fmsk); 
			isperiodic(idx,jdx)
			%p_periodic(idx,jdx) = stati.pn;
		
			b(idx,jdx) = real(bi(idx,jdx)); 
		end % for idx
		if (0)
			Sr = gampdf_man(fr,p(idx,1),p(jdx,2));
			Tr = sqrt(Sr);
			bi = ifft2(Tr.*fe);
			br(idx,:) = real(bi(idx,:));
		end
		end % for jdx
end % switch mode

% store generated pattern
save(filename,'b','Scx','ScxScy','Scy_div_Scx','Scy','p','isperiodic','p_periodic');
end % if ~exist file

% display
figure(1);
clf();
b_ = b-mean(b(:));
b_ = b_/std(b_,[],'all');
b_ = normcdf(b_);
xS = x;
yS = y;

imagesc(xS,yS,(b_')>0.5);
colormap gray
daspect([1,2,1])
axis equal;
axis tight
axis xy
%Sc_tick = [0.05,0.1,0.2,0.5,1,2,5,10,20];
Sc_tick = 2.^(-4:5);
x_tick = interp1(log10(ScxScy),x,log10(Sc_tick),'linear','extrap');
set(gca,'xtick',x_tick,'xticklabel',num2str(cvec(round(Sc_tick,3))))
xlabel('Regularity$_{2\mathrm{d}}$ $\displaystyle\frac{S_{cx}\cdot S_{cy}}{\lambda_c^2}$','interpreter','latex')
switch (mode)
case {'equal'}
case {'independent'}
	Sc_tick = [0.2,0.5,1,2,5,10,20];
	y_tick = linspace(0,L(1),7);
	set(gca,'ytick',y_tick','yticklabel',num2str(cvec(round(Sc_tick,1))))
case {'relative'}
	y_tick = [0.1,0.2,0.5,1,2,5,10];
	set(gca,'ytick',log10(y_tick'),'yticklabel',num2str(cvec(round(y_tick,1))))
	ylabel('Anisotropy of 2d-regularity S_{cy}/S_{cx}');
end % switch mode

Sc_tick = 2.^(-4:5);
y_tick = interp1(log10(Scy_div_Scx),y,log10(Sc_tick),'linear','extrap');
set(gca,'ytick',(y_tick'),'yticklabel',num2str(cvec(round(Sc_tick,2))))
ylabel('Anisotropy of Regularity$_{2\mathrm{d}}$ $S_{cy}/S_{cx}$\hspace*{3em}','interpreter','latex');

hold on
dat = load('mat/patterns-metastudy.mat');
lc = cvec(1./arrayfun(@(x) x.fc.x.hp,dat.stat));
ismodel = cvec(arrayfun(@(x) x.ismodel,dat.stat));
isisotropic = cvec(arrayfun(@(x) x.isisotropic,dat.stat));
exclude = cvec(arrayfun(@(x) x.exclude,dat.stat));
Scx_ = cvec(arrayfun(@(x) x.Sc.x.hp,dat.stat));
Scy_ = cvec(arrayfun(@(x) x.Sc.y.hp,dat.stat));

Scx_Scy_     = Scx_.*Scy_;
Scy_div_Scx_ = Scy_./Scx_;

col = 'br';
for idx=1:2
	fdx = ismodel == (idx-1) & isisotropic == 0 & exclude == 0;
	qxqy = quantile(cvec(Scx_Scy_(fdx))./(lc(fdx).^2),[0.25,0.5,0.75]); 
	qy_d_qx = quantile(cvec(Scy_div_Scx_(fdx)),[0.25,0.5,0.75]); 
	x_ = interp1(log(ScxScy),x,log10(qxqy),'linear');
	y_ = interp1(log(Scy_div_Scx),y,log10(qy_d_qx),'linear');
	errorbar(x_(2),y_(2),y_(2)-y_(1),y_(3)-y_(1),x_(2)-x_(1),x_(3)-x_(2),[col(idx)],'linewidth',2)
end % for idx

colormap(flipud(colormap_vegetation()))

figure(10);
imagesc(xS,yS,p_periodic');
%isperiodic');
colormap gray
daspect([1,2,1])
axis equal;
axis tight;
axis xy;
set(gca,'xtick',x_tick,'xticklabel',num2str(cvec(round(Sc_tick,3))))
xlabel('Regularity$_{2\mathrm{d}}$ $\displaystyle\frac{S_{cx}\cdot S_{cy}}{\lambda_c^2}$','interpreter','latex')
set(gca,'ytick',(y_tick'),'yticklabel',num2str(cvec(round(Sc_tick,2))))         

if (pflag)
	pdfprint(1,'img/pattern-synthetic-anisotropic-sweep.pdf',ps)
end

