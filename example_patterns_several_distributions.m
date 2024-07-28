% Mon  1 Jul 12:14:00 CEST 2024
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
% plot three patterns with same regularity but different distributions
% for illustration
if (~exist('pflag','var'))
	pflag =0;
end
fflag = pflag;

% characteristic frequency
fc = 1;
% regularity
reg = 2;
% density maximum
Sxpc = reg / fc;
Syc = Sxpc;
% spatial extent
L  = 10/fc;
% spatial resolution
dx = fc/10;
% spectral resolution
dfx = 1./L;
dfy = 1./L;
% number of grid points
n = L/dx;

% spatial axis
x = (0:n-1)*L/n;
y = x;
% spectral axis
fx = fourier_axis(x);
fy = fourier_axis(y);
fx = cvec(fx);
fy = rvec(fy);

distribution_C = {'gauss', 'laplace', 'cauchy'};

% for each distribution
for ddx=1:length(distribution_C)
distribution = distribution_C{ddx};

rng(0)
reg = [];
reg_ = [];
	% generate random pattern
	e  = randn(n);
	switch (distribution)
	case {'gauss'}
		% density parameter
		[f0x,sx] = normalwrappedpdf_mode2par(fc,Sxpc);
		[f0y,sy] = normpdf_mode2par(0,Syc);
		% density along axis
		Sx_ = 0.5*normalwrappedpdf(fx,f0x,sx);
		Sy_ = normpdf(fy,f0y,sy);
		% two-dimensional density
		Sxy_ = cvec(Sx_)*rvec(Sy_);
	case {'laplace'}
		[f0x,sx] = laplacewrappedpdf_mode2par(fc,Sxpc);
		[f0y,sy] = laplacepdf_mode2par(0,Syc);
		Sx_ = 0.5*laplacewrappedpdf(fx,f0x,sx);
		Sy_ = laplacepdf(fy,f0y,sy);
		Sxy_ = cvec(Sx_)*rvec(Sy_);
		sc = 0.785*0.8;
		sx = sc*sx;
		sy = sc*sy;
		fr1   = hypot((cvec(fx)-f0x)/sx,rvec(fy)/sy);
		fr2   = hypot((cvec(fx)+f0x)/sx,rvec(fy)/sy);
		Sxy_     = 1/(sx*sy)*(exp(-fr1) + exp(-fr2));
	case {'cauchy'}
		[fx0,sx] = cauchywrappedpdf_mode2par(fc,Sxpc);
		[fy0,sy] = cauchypdf_mode2par(0,Syc);
		Sx_      = 0.5*cauchywrappedpdf(fx,fx0,sx);
		Sy_      = cauchypdf(fy,f0y,sy);
		% note, this misses the normalization
		sc = 0.25; % for Sc = 1, 0.2518
		sc = 0.245; % for Sc = 1.4
		sx = sc*sx;
		sy = sc*sy;
		Sxy_       = 0.5*( 1./(1 + ((fx-fx0)./sx).^2 + ((fy-fy0)./sy).^2) ...
			   + 1./(1 + ((fx+fx0)./sx).^2 + ((fy-fy0)./sy).^2) );
		% S_ = cvec(Sx_)*rvec(Sy_);
	end
		%S_ = 2*S_/(sum(S_,'all')*dfy*dfx);
		% normalize
		fdx = fx>=0;
		Sxy_ = Sxy_/(sum(Sxy_(fdx,:),'all')*dfy*dfx);
		Sx = sum(Sxy_,2)*dfy;
		Sy = sum(Sxy_,1)*dfx;
[2*max(Sx_) max(Sy_) max(Sxy_,[],'all')]

% autocorrelation
R = real(ifft2(Sxy_));
% pattern
b = ifft2(sqrt(Sxy_).*fft2(e));
b = real(b);
sd = std(b,[],'all');


figure(1);
subplot(2,3,ddx);
plot(fx,Sx)
axis square
axis tight

subplot(2,3,ddx+3);
plot(fy,Sy);
axis square
axis tight

splitfigure([4,3],[1,ddx],fflag);
imagesc(x,y,b')
caxis(quantile(b(:),[0.625,0.6875]));
axis square
axis tight
xlabel('Distance x/\lambda_c');
ylabel('Distance y/\lambda_c');
colormap(colormap_vegetation(256));
title([upper(distribution(1)),distribution(2:end)]);

splitfigure([4,3],[1,3+ddx],fflag);
% imagesc(t);
b_ = medfilt1(b,round(1/dx),[],2);
t=b_>quantile(b(:),0.5);
m = round(0.5/dx);
t=t(:,m+1:end-m);
% rising limbs
t = t(2:end,:) & (t(2:end,:) ~= t(1:end-1,:));
[id,jd] = find(t);
l=diff(id);
l=l(l>0);
histogram(dx*l(:),(1:3./dx)*dx,'normalization','pdf')

if (0)
%  clf; b=b>quantile(b(:),0.7); st=strel('disk',1); b_=imopen(b,st); subplot(2,2,1); imagesc(b); subplot(2,2,2); imagesc(b_); axis equal;

splitfigure([4,3],[1,3+ddx],fflag);
q = 0.75;
%b_ = b>quantile(b,q,'all');
%imagesc(x,y,b_')
rw = 2.5;
xw = (-3:3);
w  = hypot(cvec(xw),rvec(xw))<=rw
b = ordfilt2(b,round(sum(w(:))/2),w);
%b = medfilt2(b,[51,51]);
imagesc(x,y,b')
caxis(quantile(b(:),[0.625,0.6875]));

colormap(colormap_vegetation(256));
axis equal
axis tight
xlabel('Distance x/\lambda_c');
ylabel('Distance y/\lambda_c');
title([upper(distribution(1)),distribution(2:end)]);
end % if 0

% plot two dimensiona density
splitfigure([4,3],[1,6+ddx],fflag);
imagesc(fftshift(fx),fftshift(fy),fftshift(Sxy_))
axis equal
axis tight
xlim([-3,3])
ylim([-3,3])

% plot two dimensional acf
splitfigure([4,3],[1,9+ddx],fflag);
imagesc(x,y,fftshift(R));
%fftshift(fx),fftshift(fy),fftshift(S_))
axis equal
axis tight
end % for ddx

if (pflag)
	ps = 3.5;
	pdfprint(11,'img/stochastic-pattern-gauss.pdf',ps);
	pdfprint(12,'img/stochastic-pattern-laplace.pdf',ps);
	pdfprint(13,'img/stochastic-pattern-cauchy.pdf',ps);
end
