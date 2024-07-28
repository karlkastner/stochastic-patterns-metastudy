% Sat 27 Jul 20:08:29 CEST 2024
% TODO explore different methods for fitting

fc = 1;
Lx = 10/fc;
Ly = 10/fc;
dx = fc/10;
dy = fc/10;
Sxpc = logspace(-1,1,10);
Syc = Sxpc;
nx = Lx/dx;
ny = Ly/dy;
dfx = 1./Lx;
dfy = 1./Ly;

fx = fourier_axis(Lx,nx); 
fy = fourier_axis(Ly,ny); 

dist_C = {'normal','lognormal','cauchy','laplace'}

% TODO sanity check of mode computation 
rng(0)

% number of mote carlo reps
m = 1000;
if (~exist('reg','var'))
reg_ = [];
leg_C = {'direct','normal','lognormal','cauchy','laplace'}
for idx=1:length(Sxpc)
    for jdx=1:length(dist_C)
	switch (dist_C{jdx})
	case {'normal'}
	[a,b] = normalmirroredpdf_mode2par(fc,0.5*Sxpc(idx));
	[fc__(idx,1), Sxc__] = normalmirroredpdf_mode(a,b);
	Sxpc__(idx,1) = 2*Sxc__;
	Sx = normalmirroredpdf(fx,a,b);
	[a,b] = normpdf_mode2par(fc,Syc(idx));
	Sy = normpdf(fy,a,b);
	Sxy = cvec(Sx)*rvec(Sy);
	

	case {'lognormal'}
	[a,b] = lognmirroredpdf_mode2par(fc,0.5*Sxpc(idx));
	[fc__(idx,2), Sxpc__(idx,2)] = lognpdf_mode(a,b);
	Sx    = lognmirroredpdf(fx,a,b);
	[a,b] = normpdf_mode2par(fc,Syc(idx));
	Sy    = normpdf(fy,a,b);
	Sxy   = cvec(Sx)*rvec(Sy);
	case {'cauchy'}
	[a,b] = cauchymirroredpdf_mode2par(fc,0.5*Sxpc(idx));
	[fc__(idx,3), Sxc___] = cauchymirroredpdf_mode(a,b);
	Sxpc__(idx,3) = 2*Sxc__;
	Sx    = cauchymirroredpdf(fx,a,b);
	[a,b] = cauchypdf_mode2par(fc,Syc(idx));
	Sy    = cauchypdf(fy,a,b);
	Sxy = cvec(Sx)*rvec(Sy);
	case {'laplace'}
	[a,b] = laplacemirroredpdf_mode2par(fc,0.5*Sxpc(idx));
	[a,b] = laplacemirroredpdf_mode2par(fc,0.5*Sxpc(idx));
	[fc__(idx,4), Sxc___] = laplacemirroredpdf_mode(a,b);
	Sxpc__(idx,4) = 2*Sxc__;
	Sx    = laplacemirroredpdf(fx,a,b);
	[a,b] = laplacepdf_mode2par(fc,Syc(idx));
	Sy    = laplacepdf(fy,a,b);
	Sxy   = cvec(Sx)*rvec(Sy);
	end
	T = sqrt(Sxy);
        for kdx=1:m
		% generate pattern
		e = randn(nx,ny);
		b = ifft2(T.*fft2(e));
		hatSxy = abs(fft2(b)).^2;
		hatSxy = hatSxy/(sum(hatSxy,'all')*dfx*dfy);

		barSx = sum(hatSxy,2).*dfy;
		% estimate density conventionally
		[Sxc_,mdx] = max(barSx);
		fc_ = abs(fx(mdx));
		reg(idx,jdx,1,kdx) = 2*Sxc_*fc_;
		[Sxc_, fc_] = extreme3(fx,barSx,mdx);
		reg(idx,jdx,6,kdx) = 2*Sxc_*fc_;

		% fit mirrored-normal
		[par0(1),par0(2)] = normalmirroredpdf_mode2par(fc,0.5*Sxpc(idx));
		par = lsqnonlin(@(par) sqrt(barSx) - sqrt(normalmirroredpdf(fx,par(1),par(2))),par0);
		[fc_, Sxc_] = normalmirroredpdf_mode(par(1),par(2));
		reg(idx,jdx,2,kdx) = 2*Sxc_*fc_;

		% fit log-normal
		[par0(1),par0(2)] = lognmirroredpdf_mode2par(fc,0.5*Sxpc(idx));
		par = lsqnonlin(@(par) sqrt(barSx) - sqrt(lognmirroredpdf(fx,par(1),par(2))),par0);
		%[fc_,Sc_] = lognmirroredpdf_mode(par(1),par(2));
		[fc_,Sxpc_] = lognpdf_mode(par(1),par(2));
		reg(idx,jdx,3,kdx) = Sxpc_*fc_;
		% fit cauchy lorentz
		% fit laplace
		% estimate reg from par

		% fit cauchy
		[par0(1),par0(2)] = cauchymirroredpdf_mode2par(fc,0.5*Sxpc(idx));
		par = lsqnonlin(@(par) sqrt(barSx) - sqrt(cauchymirroredpdf(fx,par(1),par(2))),par0);
		[fc_,Sc_] = cauchymirroredpdf_mode(par(1),par(2));
		reg(idx,jdx,4,kdx) = 2*Sc_*fc_;

		% fit laplace
		% TODO, we should account for the tail that is truncated
		[par0(1),par0(2)] = laplacemirroredpdf_mode2par(fc,0.5*Sxpc(idx));
		par = lsqnonlin(@(par) sqrt(barSx) - sqrt(laplacemirroredpdf(fx,par(1),par(2))),par0);
		[Sc_,fc_] = laplacemirroredpdf_mode(par(1),par(2));
		reg(idx,jdx,5,kdx) = 2*Sc_/fc_;
  end
 end
end

res = reg - cvec(Sxpc);
rmse = rms(res,4);
end
isgood = double(abs(fc__ - fc) < 1e-2);
isgood(~isgood) = NaN;
rmse_ = rmse.*isgood;

figure(1)
for idx=1:length(dist_C)
	subplot(2,2,idx);
	loglog(Sxpc,squeeze(rmse_(:,idx,:)));
	title(dist_C{idx})
	if (1==idx)
		legend(leg_C);
	end
	xlim(limits(Sxpc*fc));
end
figure(2)
for idx=1:length(dist_C)
	subplot(2,2,idx);
	loglog(Sxpc,squeeze(rmse_(:,idx,:))./rmse_(:,idx,1));
	title(dist_C{idx})
	if (1==idx)
		legend(leg_C);
	end
	xlim(limits(Sxpc*fc));
end
Sxpc__
fc__



