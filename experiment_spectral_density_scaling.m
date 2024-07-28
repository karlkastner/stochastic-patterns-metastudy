% 2022-03-30 13:29:31.635728730 +0200
if (~exist('pflag','var'))
	pflag = 0;
end
fflag = pflag;
% plotscale
ps = 3.5;
% linewidth
lw = 1.5;
d=0;

% spatial extent
L = 100; 
% spectral resolution
dfx = 1/L;
% spatial axis
x=linspace(0,L,1e4)';
% spectral axis
fx = fourier_axis(x);

for jdx=1

%if (1==jdx)
%	% regularity
%	reg = [0.5,1,1.5,1,1];
%	% characteristic length scale
%	lc  = [1,1,1,0.5,2];
%	% density maxima
%	Sc  = reg.*lc;
%else
%	reg = [0.5,1,2];
%	lc = [1,1,1];
%	Sc = reg.*lc;
%end
%

% characteristics length scale
lc  = [0.5,0.5,0.5,2,2,2];
% characteristic frequncy
fc = 1./lc;
% regularity
reg = [0.5,1,2,0.5,1,2];
% density maxima
Sxpc = reg.*lc;

Sxp = [];
for idx=1:length(reg)
	[p,q] = lognpdf_mode2par(fc(idx),Sxpc(idx));
	Sxp(:,idx) = lognpdf(fx,p,q);
	leg_C{idx} = sprintf('\\lambda_c=%0.1f S_c=%1.1f',lc(idx),Sxpc(idx))
	%leg_C{idx} = sprintf('\\lambda_c=%0.1f S_c=%1.1f Reg=%1.1f',lc(idx),Sc(idx),reg(idx))
end
% normalize
Sxp = Sxp ./ (sum(Sxp)*dfx);
fdx=fx>0;

splitfigure([2,2],[1,1+2*(jdx-1)],fflag);
cla();
%ls_C = {'k--','k-','k-.','r-','b-','-'}
ls_C = {'b-','b-.','b--','r-','r-.','r--'}
fdx = fx>0;
for idx=1:length(reg)
	plot(fx(fdx),Sxp(fdx,idx),ls_C{idx},'linewidth',lw);
	hold on
end
%plot(fx(fdx),S(fdx,2),'r','linewidth',lw);
%plot(fx(fdx),S(fdx,3),'k','linewidth',lw);
%xlabel('Wavenumber $$k \, \frac{[L]}{2 \pi}$$','interpreter','latex')
xlabel('Wavenumber $$k \, [L] / (2 \pi)$$','interpreter','latex')
ylabel('Density $S_x^+ / [L]$','interpreter','latex','rot',90);
xlim([0,4])
legend(leg_C{:})

%subplot(2,2,2+2*(jdx-1));
splitfigure([2,2],[1,2+2*(jdx-1)],fflag);
cla
%ls_C = {'b-','r--','k:'}
ls_C = {'k-','k-.','k--'}; %,'r--','r-','r-.'}
for idx=1:3
plot(fx(fdx)/fc(idx),Sxp(fdx,idx)*fc(idx)-d,ls_C{idx},'linewidth',lw);
leg_C{idx} = sprintf('S_c/\\lambda_c=%0.1f',reg(idx))
hold on
end
%plot(fx(fdx)/fc,S(fdx,1)*fc,'r--','linewidth',lw);
%plot(fx(fdx)/fc,S(fdx,3)*fc,'k:','linewidth',lw);
xlabel('Wavenumber $k / k_c$','interpreter','latex')
ylabel('Density $S_{x}^+ / \lambda_c$','interpreter','latex','rot',90);
xlim([0,4])
ylim([0,4])
legend(leg_C)
%title('Normalized');
end

if (pflag)
	pdfprint(11,'img/regularity-parameter.pdf',ps);
	pdfprint(12,'img/regularity-parameter-normalized.pdf',ps);
end

