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
L = 40;
m = L;
n = m*L;
% s=[2.8,2,1.4,1,0.7,0.5];
rng_ = 1;
filename=sprintf('img/synthetic-isotropic-pattern-rng-%d-L-%d-m-%d.mat',rng_,L,m);
if (exist(filename,'file'))
	load(filename);
else
s   = flipud(2.^(-2:0.5:3.5)');
bb  = [];
Scr = [];
Sca = [];
for idx=1:length(s)
 disp(idx)
 rng(rng_);
%if (0)
 [b,x,y]=hexagonal_pattern(1,n,L,0,0,s(idx),0.5,1);
%else
%	b =bb(:,:,idx);
%end
 bb(:,:,idx) = b;
 S=abs(fft2(b-mean(b(:)))).^2;
 Sbar=trifilt2(S,3);
 [Sr,fr] = periodogram_radial(S,[L,L]);
 fdx=find(fr==1);
 Scr(idx,1) = Sr.normalized(fdx);
 Scr(idx,2) = max(Sr.normalized)
 [Sa,angle] = periodogram_angular(S,[L,L]);
 Sca(idx,1) = interp1(angle,Sa,0,'linear');
 Sca(idx,2) = max(Sa);
figure(4);
subplot(4,6,idx)
plot(angle,Sa);
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
end
end;

%caxis(quantile(b(:),[0.4,0.6]));
disp('interpolating');
Sc_ = sort(Scr(:,1));
%iSc = logspace(log10(Sc_(1)),log10(Sc_(end)),n);
%iSc = logspace(-1,1,n);
iSc = logspace(log10(0.125),log10(16),n);
bi = [];
for idx=1:n
 for jdx=1:n
  bi(idx,jdx) = interp1(Sc_,squeeze(bb(idx,jdx,:)),iSc(jdx),'linear');
 end
end
%if (0)
save(filename,'x','y','Scr','Sca','s','bb','bi');
end
figure(6)
 clf;plot(1./s,[Scr,Sca]) 


%end
%figure(4); imagesc(bi>0.5); axis equal; axis tight
figure(5);
 imagesc(x,y,bi>0.6);
 axis equal;
 axis tight;
 Sct = 2.^(-3:4);
 xt=interp1(log2(iSc),x,log2(Sct),'linear');
 set(gca,'xtick',xt,'xticklabel',num2str(cvec(Sct)))
colormap(flipud(gray));
xlabel('Regularity S_{cr}/\lambda_c');
ylabel('Distance y/\lambda_c');

Lp = 10:10:L;
for idx=1:length(Lp)
	figure(5);
	ylim(Lp(idx)/2*[-1,1]);
	pdfprint(5,sprintf('img/synthetic-isotropic-pattern-L-%d-m-%d-Lp-%d.pdf',L,m,Lp(idx)),2);
end

if (0)

rng(1);
% 2
 s=0.003;
s = 0.001;
n=1001;
 L=40;
 x = linspace(-L,L,n)';
 y = x';
 z0 = cos(2*pi*x+0.*y);
 z=0;
 fx=fourier_axis(x);
 fy=fx';
 %Sx=spectral_density_bandpass_continuous(fx,1,200);
 %Sy = exppdf(abs(fy),0.001);
 %S=Sx*Sy;
if (1)
 [fx,fy,fr,a] = fourier_axis_2d([L,L],[n,n]);
 fy = rvec(fy);
% S=(normpdf(fx,1,s)+normpdf(fx,-1,s)).*normpdf(fy,0,4*s);
 %S = normpdf(fr,1,s).*(normpdf(a,0,s)+normpdf(a,pi,s)+normpdf(a,-pi,s));

d = hypot(fx-1,fy);
d2 = hypot(fx+1,fy);
S = exp(-d/sqrt(s)) + exp(-d2/sqrt(s));

 T=sqrt(S);
% e=randn(n);
 for idx=0:2;
 e=randn(n);
 f=fft2(e);
% f=abs(f);
 z0 = real(ifft2(T.*f));
  z=z+imrotate(z0,60*idx,'bilinear','crop');
 end;
else
[fx,fy,fr,a] = fourier_axis_2d([L,L],[n,n]);
 fy = rvec(fy);
S0 = normpdf(fr,1,s).*(normpdf(a,0,s)+normpdf(a,pi,s)+normpdf(a,-pi,s));
S0=(normpdf(fx,1,s)+normpdf(fx,-1,s)).*normpdf(fy,0,4*s);

d = hypot(fx-1,fy);
d2 = hypot(fx+1,fy);
S0 = exp(-d/sqrt(s)) + exp(-d2/sqrt(s));

S = 0;
for idx=0:2
 S = S+fft_rotate(S0,idx*60);
end
 T=sqrt(S);
 e=randn(n);
 f=fft2(e);
 %f=dct2(e);
 %f=abs(f);
 Shat = (T.*f);
 z = real(ifft2(Shat));
 %z = real(idct2(Shat));
end
 z(hypot(x,y)>L)=NaN;
 imagesc(z/nanstd(z,[],'all')>0.75);
%imagesc(z)
 axis equal

end
