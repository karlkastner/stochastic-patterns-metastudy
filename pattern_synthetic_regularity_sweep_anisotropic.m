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
pdfx_str  = 'normal';
%pdfx_str  = 'lognormal';
%pdfy_str = 'laplace';
pdfy_str  = 'normal'
ps = 1.5;
pc = 0.375;
cmap = colormap_vegetation2();
%colormap(flipud((1-pc)*gray() + pc*colormap_vegetation()))
% mode='rely';
% TODO laplace, lorentzian

% distribution along the primary axis
switch (pdfx_str)
case {'normal'}
	pdfx          = @normalwrappedpdf;
	pdfx_mode2par = @(fc,Sc) normalwrappedpdf_mode2par(fc,0.5*Sc);
case {'lognormal'}
	pdfx          = @longmirroredpdf;
	pdfx_mode2par = @lognpdf_mode2par;
end

% distribution along the secondary axis
switch (pdfy_str)
case {'normal'}
	pdfy = @normpdf;
	pdfy_mode2par = @normpdf_mode2par;
case {'laplace'}
	pdfy = @laplacepdf;
	pdfy_mode2par = @laplacepdf_mode2par;
end

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
	Sxpc = 2*logspace(-1,1,n(1))';
	Syc = logspace(-1,1,n(1));
case {'independent'}
	Sxpc = 2*logspace(-1,1,n(1))';
	L(2) = L(1);
	n(2) = n(1);
	Syc = logspace(-1,1,n(2));
case {'rely'}
	L(2) = 30;
	n(2) = m*L(2);
	Syc_rel = logspace(-1,1,n(2));
	Syc = cvec(Sxpc)*Syc_rel;
case {'relxy'}
	SxpcSyc_lim       = [0.125,16];
	Syc_div_Sxpc_lim  = [2.^-1.5,2.^1.5];
	L = [60,30];
	L(2) = L(1)*range(log10(Syc_div_Sxpc_lim))/range(log10(SxpcSyc_lim));
	m = 10;
	n  = round(m*L);
	SxpcSyc      = logspace(log10(SxpcSyc_lim(1)),log10(SxpcSyc_lim(2)),n(1))';
	Syc_div_Sxpc = logspace(log10(Syc_div_Sxpc_lim(1)),log10(Syc_div_Sxpc_lim(2)),n(2))';
	Sxpc = sqrt(SxpcSyc./rvec(Syc_div_Sxpc));
	Syc = sqrt(SxpcSyc.*rvec(Syc_div_Sxpc));

end

% output file
filebase_str = sprintf('synthetic-pattern-anisotropic-Sx-%s-Sy-%s-Sxpc-%f-%f-Syc-%f-%f-L-%d-%d-dx-%d-rng-%d',pdfx_str,pdfy_str,min(Sxpc(:)),max(Sxpc(:)),min(Syc(:)),max(Syc(:)),L(1),L(2),dx(1),rng_);
filename_str = ['mat/',filebase_str,'.mat'];
if (exist(filename_str,'file'))
	disp('Loading file');
	load(filename_str);
else

% axes in real space
x = innerspace(0,L(1),n(1));
y = innerspace(0,L(2),n(2));

% axes in frequency space
[fx,fy,fr] = fourier_axis_2d(L,n);

e  = randn(n);
fe = fft2(e);

% vary y-regularity on y-axis
%-> vary Syc on y-axis

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

%	if (1)
		[p.x(idx,jdx,1),p.x(idx,jdx,2)] = pdfx_mode2par(fc,Sxpc(idx,jdx)); 
		%if (~strcmp(pdfy_str,'exp'))
		[p.y(idx,jdx,1),p.y(idx,jdx,2)] = pdfy_mode2par(0,Syc(idx,jdx)); 
		%else
		%	[p.y(idx,jdx,1)] = pdfy_mode2par(Syc(idx,jdx)); 
		%end
%	else
%		[p.x(idx,jdx,1),p.x(idx,jdx,2)] = pdf_mode2par(fc,Sxpc(idx),[p.x(idx-1,jdx,1),p.x(idx-1,jdx,2)],flag);
%		[p.y(idx,jdx,1),p.y(idx,jdx,2)] = pdf_mode2par(0,Syc(idx),[p.y(idx-1,jdx,1),p.y(idx-1,jdx,2)],flag);
%	end
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
			Sx = pdfx(fx,p.x(idx,jdx,1),p.x(idx,jdx,2));
			%if (~strcmp(pdfy_str,'exp'))
			Sy = pdfy(fy,p.y(idx,jdx,1),p.y(idx,jdx,2));
			%else
			%	Sy = pdfy(abs(fy),p.y(idx,jdx,1));
			%end
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
			%isperiodic(idx,jdx)
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
save(filename_str,'x','y','L','dx','b','Sxpc','SxpcSyc','Syc_div_Sxpc','Syc','p','isperiodic','p_periodic');
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
x_tick = interp1(log10(SxpcSyc),x,log10(Sc_tick),'linear','extrap');
set(gca,'xtick',x_tick,'xticklabel',num2str(cvec(round(Sc_tick,3))))
xlabel('Regularity$_{xy}$ $\displaystyle\frac{S_{xc}^+\cdot S_{yc}}{\lambda_c^2}$','interpreter','latex')
switch (mode)
case {'equal'}
case {'independent'}
	Sc_tick = [0.2,0.5,1,2,5,10,20];
	y_tick = linspace(0,L(1),7);
	set(gca,'ytick',y_tick','yticklabel',num2str(cvec(round(Sc_tick,1))))
case {'relative'}
	y_tick = [0.1,0.2,0.5,1,2,5,10];
	set(gca,'ytick',log10(y_tick'),'yticklabel',num2str(cvec(round(y_tick,1))))
	ylabel('Anisotropy of regularity S_{yc}/S_{xc}^+');
end % switch mode

Sc_tick = 2.^(-4:5);
y_tick = interp1(log10(Syc_div_Sxpc),y,log10(Sc_tick),'linear','extrap');
set(gca,'ytick',(y_tick'),'yticklabel',num2str(cvec(round(Sc_tick,2))))
ylabel('Anisotropy of Regularity $S_{yc}/S_{xc}^+$\hspace*{3em}','interpreter','latex');

hold on
dat = load('mat/patterns-metastudy.mat');
lc = cvec(1./arrayfun(@(x) x.fc.x.hp,dat.stat));
ismodel = cvec(arrayfun(@(x) x.ismodel,dat.stat));
isisotropic = cvec(arrayfun(@(x) x.isisotropic,dat.stat));
exclude = cvec(arrayfun(@(x) x.exclude,dat.stat));
Sxpc_ = cvec(arrayfun(@(x) x.Sc.xp.hp,dat.stat));
Syc_ = cvec(arrayfun(@(x) x.Sc.y.hp,dat.stat));

Sxpc_Syc_     = Sxpc_.*Syc_;
Syc_div_Sxpc_ = Syc_./Sxpc_;

col = 'br';
for idx=1:2
	fdx = ismodel == (idx-1) & isisotropic == 0 & exclude == 0;
	qxqy = quantile(cvec(Sxpc_Syc_(fdx))./(lc(fdx).^2),[0.25,0.5,0.75]);
	qy_d_qx = quantile(cvec(Syc_div_Sxpc_(fdx)),[0.25,0.5,0.75]); 
	x_ = interp1(log(SxpcSyc),x,log10(qxqy),'linear');
	y_ = interp1(log(Syc_div_Sxpc),y,log10(qy_d_qx),'linear');
	errorbar(x_(2),y_(2),y_(2)-y_(1),y_(3)-y_(1),x_(2)-x_(1),x_(3)-x_(2),[col(idx)],'linewidth',2)
end % for idx

colormap(cmap)

figure(10);
imagesc(xS,yS,p_periodic');
%isperiodic');
colormap gray
daspect([1,2,1])
axis equal;
axis tight;
axis xy;
set(gca,'xtick',x_tick,'xticklabel',num2str(cvec(round(Sc_tick,3))))
xlabel('Regularity$_{2\mathrm{d}}$ $\displaystyle\frac{S_{xc}\cdot S_{yc}}{\lambda_c^2}$','interpreter','latex')
set(gca,'ytick',(y_tick'),'yticklabel',num2str(cvec(round(Sc_tick,2))))         

if (pflag)
	pdfprint(1,'img/pattern-synthetic-anisotropic-sweep.pdf',ps)
	pdfprint(1,['img/',filebase_str,'.pdf'],ps)
end

