clear all
close all

rep='../RESULTS/';
filep='parametres_fault_rns.cfg';
par=readparametres_cfg(rep,filep);

kfig=0;
%--Plot spatial average variables time series-------------%
%[kfig]=plt_av_var(par,rep,kfig);
%--Plot local variables time series-----------------------%
%[kfig]=plt_loc_var(par,rep,kfig);
%--Plot HV profils----------------------------------------%
%[kfig]=plt_hv_profils(par,rep,kfig);
%--Plot maps----------------------------------------------%
[kfig]=plt_maps(par,rep,kfig);
%--Plot catalog-------------------------------------------%
%[kfig]=plt_catalog(par,rep,kfig);

function [kfig]=plt_av_var(par,rep,kfig)
%-------------------------------------------------------------%
%--Plot spatial average of variables (qmoync.m)---------------%
%-------------------------------------------------------------%


year=365*24*3600;


[t,dt,v,vmax,vmvw,th,tau,u]=extractimeseriesnc(rep);

%--compute moment rate
mrate=v*par.nx*par.dx*par.Y/(2*(1+par.nu));
mratevw=vmvw*par.nx*par.dx*par.Y/(2*(1+par.nu));

%--Plot
figure(kfig+1);
subplot(2,2,1,'align');
semilogy(t/year,v,'-k');
hold on
semilogy(t/year,vmax,'-r');
xlabel('t (years)');
ylabel('v (m/s)');

subplot(2,2,2,'align');
semilogy(t/year,th,'-k');
xlabel('t (years)');
ylabel('$\theta$','Interpreter','latex');

subplot(2,2,3,'align');
plot(t/year,u,'-k');
hold on
plot(t/year,par.vp*t,'--r');
xlabel('t (years)');
ylabel('u (m)');

subplot(2,2,4,'align');
plot(t/year,tau*1e-6,'-k');
xlabel('t (years)');
ylabel('$\tau$ (MPa)','Interpreter','latex');

figure(kfig+2);
semilogy(t/year,mrate,'-k');
hold on
semilogy(t/year,mratevw,'-r')
xlabel('t (years)','Interpreter','latex');
ylabel('Moment rate (N.m/s)','Interpreter','latex');
legend('VW+VS','VW')

figure(kfig+3)
semilogy([1:length(t)],dt,'o-k')
xlabel('t (years)','Interpreter','latex');
ylabel('Moment rate (N.m/s)','Interpreter','latex');
xlabel('t (years)','Interpreter','latex');
ylabel('time step dt (s)','Interpreter','latex');

kfig=kfig+3;

end



function [kfig]=plt_loc_var(par,rep,kfig)
%------------------------------------------------%
%--Plot local variables (qlocnc.m)---------------%
%------------------------------------------------%


year=365*24*3600;


[t,vloc,thloc,tauloc,uloc,ploc,qdarcyloc]=extract_local_timeseries_nc(rep,par);

for i=1:par.nploc

    figure(kfig+1);
    subplot(3,2,1,'align');
    semilogy(t/year,vloc(:,i),'-');
    hold on
    xlabel('t (years)','Interpreter','latex');
    ylabel('v (m/s)','Interpreter','latex');

    subplot(3,2,2,'align');
    semilogy(t/year,thloc(:,i),'-');
    hold on
    xlabel('t (years)');
    ylabel('$\theta$','Interpreter','latex');

    subplot(3,2,3,'align');
    plot(t/year,uloc(:,i),'-');
    hold on
    plot(t/year,par.vp*t,'--r');
    xlabel('t (years)','Interpreter','latex');
    ylabel('u (m)','Interpreter','latex');

    subplot(3,2,4,'align');
    plot(t/year,tauloc(:,i)*1e-6,'-');
    hold on
    xlabel('t (years)','Interpreter','latex');
    ylabel('$\tau$ (MPa)','Interpreter','latex');

    subplot(3,2,5,'align');
    plot(t/year,ploc(:,i)*1e-6,'-');
    hold on
    xlabel('t (years)','Interpreter','latex');
    ylabel('$P$ (MPa)','Interpreter','latex');

    subplot(3,2,6,'align');
    plot(t/year,qdarcyloc(:,i),'-');
    hold on
    xlabel('t (years)','Interpreter','latex');
    ylabel('$q_{Darcy}$ (m/s)','Interpreter','latex');

    %pause
end

kfig=kfig+1;
end

function [kfig]=plt_hv_profils(par,rep,kfig)
%------------------------------------------------%
%--Plot HV profiles of results (hvcprofilsnc.m)--%
%------------------------------------------------%

year=365*24*3600;
pltprof=1;
kplt=1;


%--hydraulic diffusivity
par.dhy=par.k/(par.phi*par.cf*par.etaf)
par.qinj
par.xinj
par.yinj
par.toff
%

x=[-0.5*(par.nx-1)*par.dx:par.dx:0.5*(par.nx-1)*par.dx]-0.5*par.dx;
y=[-0.5*(par.ny-1)*par.dy:par.dy:0.5*(par.ny-1)*par.dy]-0.5*par.dy;

[status,result1]=system(['ls -lrth ',rep,'/maps*.nc | wc -l']);
ns(1)=str2num(result1);


np=min(ns);

[t,v,vmax,vmvw,th,tau,u]=extractimeseriesnc(rep);

ii=[];
for i=1:np
    if vmax(i)<1e-4
        if mod(i,100)==0
            ii=[ii i];
        end
    else
        if mod(i,10)==0
            ii=[ii i];
        end
    end
end
%ii=[1:np];



nnp=length(ii);

iprofil=0;

for k=1:nnp

    i=ii(k);

    finfo=ncinfo([rep,'profilsx',num2str(i),'.nc']);
    szcat=size(finfo.Variables);

    t(k)=ncread([rep,'profilsx',num2str(i),'.nc'],'time');

    %--Profils along strike (x)
    vx=ncread([rep,'profilsx',num2str(i),'.nc'],'slip rate');
    ux=ncread([rep,'profilsx',num2str(i),'.nc'],'slip');
    thx=ncread([rep,'profilsx',num2str(i),'.nc'],'state');
    taux=ncread([rep,'profilsx',num2str(i),'.nc'],'shear stress');

    if szcat(2)>6
        px=ncread([rep,'profilsx',num2str(i),'.nc'],'pore pressure');
    else
        px=zeros(size(vx));
    end


    vmaxx(k)=max(vx);

    %--Profils along depth (y)
    vy=ncread([rep,'profilsy',num2str(i),'.nc'],'slip rate');
    uy=ncread([rep,'profilsy',num2str(i),'.nc'],'slip');
    thy=ncread([rep,'profilsy',num2str(i),'.nc'],'state');
    tauy=ncread([rep,'profilsy',num2str(i),'.nc'],'shear stress');

    if szcat(2)>6
        py=ncread([rep,'profilsy',num2str(i),'.nc'],'pore pressure');
    else
        py=zeros(size(vy));
    end


    vmaxy(k)=max(vy);





    if pltprof==1
        if (mod(k,kplt)==0 | k==1)

            iprofil=iprofil+1;

            %--compute analytical approximate solution
            %panaly=press_solution(par.qinj,par.dhy,x,t(i),par.toff,par.cf,par.phi);

            figure(kfig+1);
            subplot(5,2,1,'align')
            if iprofil>1
                hold on
            end
            plot(x,px*1e-6,'-k')
            %hold on
            %plot(x,panaly,'--r')
            xlabel('x (m)');ylabel('P (MPa)');
            %pause

            subplot(5,2,2,'align')
            if iprofil>1
                hold on
            end
            plot(y,py*1e-6,'-k')
            %hold on
            %plot(x,panaly,'--r')
            xlabel('y (m)');ylabel('P (MPa)');

            subplot(5,2,3,'align')
            if iprofil>1
                hold on
            end
            plot(x,taux*1e-6,'-k')
            xlabel('x (m)');ylabel('tau (MPa)');

            subplot(5,2,4,'align')
            if iprofil>1
                hold on
            end
            plot(y,tauy*1e-6,'-k')
            xlabel('y (m)');ylabel('tau (MPa)');

            subplot(5,2,5,'align')
            if iprofil>1
                hold on
            end
            semilogy(x,vx,'-k')
            xlabel('x (m)');ylabel('v (m/s)');

            subplot(5,2,6,'align')
            if iprofil>1
                hold on
            end
            semilogy(y,vy,'-k')
            xlabel('y (m)');ylabel('v (m/s)');

            subplot(5,2,7,'align')
            if iprofil>1
                hold on
            end
            semilogy(x,thx,'-k')
            xlabel('x (m)');ylabel('$\theta$ (s)','Interpreter','latex');

            subplot(5,2,8,'align')
            if iprofil>1
                hold on
            end
            semilogy(y,thy,'-k')
            xlabel('y (m)');ylabel('$\theta$ (s)','Interpreter','latex');

            subplot(5,2,9,'align')
            if iprofil>1
                hold on
            end
            plot(x,ux,'-k')
            xlabel('x (m)');ylabel('u (m)');

            subplot(5,2,10,'align')
            if iprofil>1
                hold on
            end
            plot(y,uy,'-k')
            xlabel('y (m)');ylabel('u (m)');





        end
    end





end

kfig=kfig+1;

end

function [kfig]=plt_maps(par,rep,kfig)
%-----------------------------------------%
%--Plot 2D maps of results (mapsnc_3d.m)--%
%-----------------------------------------%

year=365*24*3600;
pltmap=1;
kplt=1;


%--hydraulic diffusivity
par.dhy=par.k/(par.phi*par.cf*par.etaf)
par.qinj
par.xinj
par.yinj
par.toff


x=[-0.5*(par.nx-1)*par.dx:par.dx:0.5*(par.nx-1)*par.dx]-0.5*par.dx;
y=[-0.5*(par.ny-1)*par.dy:par.dy:0.5*(par.ny-1)*par.dy]-0.5*par.dy;

[status,result1]=system(['ls -lrth ',rep,'/maps*.nc | wc -l']);
ns(1)=str2num(result1);


np=min(ns);




ii=[1:np];



nnp=length(ii);

mu=zeros(par.ny,par.nx);
mtau=zeros(par.ny,par.nx);
mv=zeros(par.ny,par.nx);
msigmae=zeros(par.ny,par.nx);
mqdarcy=zeros(par.ny,par.nx);
for k=1:nnp

    i=ii(k);

    finfo=ncinfo([rep,'maps',num2str(i),'.nc']);
    szcat=size(finfo.Variables);


    vx=ncread([rep,'maps',num2str(i),'.nc'],'slip rate');
    ux=ncread([rep,'maps',num2str(i),'.nc'],'slip');
    thx=ncread([rep,'maps',num2str(i),'.nc'],'state');
    taux=ncread([rep,'maps',num2str(i),'.nc'],'shear stress');

    if szcat(2)>6
        px=ncread([rep,'maps',num2str(i),'.nc'],'pore pressure');
    else
        px=zeros(size(vx));
    end

    if szcat(2)>7
        qx=ncread([rep,'maps',num2str(i),'.nc'],'darcy velocity');
    else
        qx=zeros(size(vx));
    end

    t(k)=ncread([rep,'maps',num2str(i),'.nc'],'time');
    vmax(k)=max(vx);


    mv=reshape(vx,par.ny,par.nx);
    mu=reshape(ux,par.ny,par.nx);
    mth=reshape(thx,par.ny,par.nx);
    mtau=reshape(taux,par.ny,par.nx);
    msigmae=par.smap-reshape(px,par.ny,par.nx);
    mqdarcy=reshape(qx,par.ny,par.nx);




    if pltmap==1
        if mod(k,kplt)==0

            %--compute analytical approximate solution
            %panaly=press_solution(par.qinj,par.dhy,x,t(i),par.toff,par.cf,par.phi);

            figure(1+kfig);clf
            subplot(3,2,1,'align')
            imagesc(x,y,msigmae*1e-6)
            cb=colorbar;
            xlabel('x (m)');ylabel('y (m)')
            ylabel(cb,'$\sigma-P$ (MPa)','Interpreter','latex');

            subplot(3,2,2,'align')
            imagesc(x,y,mtau*1e-6)
            cb=colorbar;
            xlabel('x (m)');ylabel('y (m)')
            ylabel(cb,'$\tau$ (MPa)','Interpreter','latex');

            subplot(3,2,3,'align')
            imagesc(x,y,mu)
            cb=colorbar;
            xlabel('x (m)');ylabel('y (m)')
            ylabel(cb,'$\delta$ (m)','Interpreter','latex');

            subplot(3,2,4,'align')
            imagesc(x,y,log10(mv))
            cb=colorbar;
            xlabel('x (m)');ylabel('y (m)')
            ylabel(cb,'log v (m/s)','Interpreter','latex');

            subplot(3,2,5,'align')
            imagesc(x,y,log10(mth))
            cb=colorbar;
            xlabel('x (m)');ylabel('y (m)')
            ylabel(cb,'log $\theta$ (m/s)','Interpreter','latex');

            subplot(3,2,6,'align')
            imagesc(x,y,mqdarcy)
            cb=colorbar;
            xlabel('x (m)');ylabel('y (m)')
            ylabel(cb,'$v_{Darcy}$ (m/s)','Interpreter','latex');







            %pause
        end
    end




end

kfig=kfig+1;

end

%%
function [kfig]=plt_catalog(par,rep,kfig)
%--------------------------------------%
%--plot catalogue (readcataloguenc.m)--%
%--------------------------------------%

year=365*24*3600;


vx=[-0.5*(par.nx-1)*par.dx:par.dx:0.5*(par.nx-1)*par.dx];
vy=[-0.5*(par.ny-1)*par.dy:par.dy:0.5*(par.ny-1)*par.dy];


%--catalogue
[t0,dt0,tev,xi,yi,xb,yb,s,m0,dtau,taui,tauf,du]=extractcatnc(rep,par);

%--magnitude
mw=(2/3)*log10(m0)-6;

%--radial distance
rad=sqrt(xb.^2 + yb.^2);

%--earthquake size
if par.dx*par.dy>0
    reqk=sqrt(s);
else
    reqk=0.5*s/par.dx;
end

%--number of earthquakes
neqk=length(t0);


figure(1+kfig);
subplot(2,2,1,'align')
plot(t0/year,mw,'ok')
xlabel('time $t$ (years)','Interpreter','latex')
ylabel('magnitude $M_w$','Interpreter','latex')

subplot(2,2,2,'align')
plot(t0/year,-dtau*1e-6,'ok')
xlabel('time $t$ (years)','Interpreter','latex')
ylabel('stress drop $\Delta \tau$ (MPa)','Interpreter','latex')

subplot(2,2,3,'align')
plot(t0/year,du,'ok')
xlabel('time $t$ (years)','Interpreter','latex')
ylabel('coseismic slip $\Delta u$ (m)','Interpreter','latex')

subplot(2,2,4,'align')
plot(t0/year,tev,'ok')
xlabel('time $t$ (years)','Interpreter','latex')
ylabel('event duration $t_e$ (s)','Interpreter','latex')

figure(2+kfig);
if par.ny>1  %--3D computation
    for i=1:neqk
        plot(xi(i),yi(i),'+k');
        hold on
        r(1)=rectangle('Position',[xb(i)-reqk(i) yb(i)-reqk(i) 2*reqk(i) 2*reqk(i)],'Curvature',1)
        set(r(1),'EdgeColor','k')
        plot(xb(i),yb(i),'+r');
        r(2)=rectangle('Position',[xb(i)-reqk(i) yb(i)-reqk(i) 2*reqk(i) 2*reqk(i)],'Curvature',1)
        set(r(2),'EdgeColor','r')
    end
    xlim([min(vx) max(vx)]);
    ylim([min(vy) max(vy)]);
    %axis equal
    xlabel('x (m)');ylabel('y (m)');

else   %--2D computation
    for i=1:neqk
        plot([xb(i)-reqk(i) xb(i)+reqk(i)],[t0(i) t0(i)]/year,'-k','Linewidth',2)
        if i==1
            hold on
        end
    end
    xlim([min(vx) max(vx)]);
    axis equal
    xlabel('x (m)');ylabel('time t (years)');
end

figure(3+kfig);
loglog(m0,tev,'ok')
xlabel('seismic moment $M_0$ (N.m)','Interpreter','latex')
ylabel('event duration $t_e$ (s)','Interpreter','latex')

kfig=kfig+3;

end




%%

%--------------------------------------%

function p=press_solution(q0,alpha,x,t,toff,beta,phi)

n=length(x);
if ((t>0) && (t<toff))
    eta=abs(x)/(2*sqrt(alpha*t));
    G=sqrt(t)*(exp(-eta.^2)/sqrt(pi)-eta.*erfc(eta));
    p=q0*G/(beta*phi*sqrt(alpha));
elseif t>toff
    eta1=abs(x)/(2*sqrt(alpha*t));
    G1=sqrt(t)*(exp(-eta1.^2)/sqrt(pi)-eta1.*erfc(eta1));
    eta2=abs(x)/(2*sqrt(alpha*(t-toff)));
    G2=sqrt(t-toff)*(exp(-eta2.^2)/sqrt(pi)-eta2.*erfc(eta2));
    p=q0*(G1-G2)/(beta*phi*sqrt(alpha));
else
    p=zeros(1,n);
end

end

%--------------------------------------%
function par=readparametres_cfg(rep,filep)

par.rep=rep;
file=[par.rep,filep];

fid=fopen(file,'r'); C=textscan(fid,'%s','Delimiter','\n','Whitespace',''); fclose(fid);
if ~startsWith(C{1}{1},'#'), C{1}{1}=['#' C{1}{1}]; end
fid=fopen(file,'w'); fprintf(fid,'%s\n',C{1}{:}); fclose(fid);

R=ini2struct(file);

par.simlab=char(R.simlab);
par.slipmode=str2num(char(R.slip_mode));
par.Y=str2num(char(R.young_mod));
par.nu=str2num(char(R.poisson_ratio));
par.rho=str2num(char(R.rho_r));
par.H=str2num(char(R.slab_thickness));
par.ess=str2num(char(R.ext_shear_stressing));
par.mu0=str2num(char(R.ref_fric_coeff));
par.alpha=str2num(char(R.dl_coeff));
par.vp=str2num(char(R.vplate));
par.v0=str2num(char(R.vstar));
par.vsis=str2num(char(R.vsis));
par.nasp=str2num(char(R.nasp));
par.aasp=str2num(char(R.a));
par.basp=str2num(char(R.b));
par.dcasp=str2num(char(R.dcasp));
par.xcasp=str2num(char(R.xcasp));
par.ycasp=str2num(char(R.ycasp));
par.rasp=str2num(char(R.rasp));
par.sigasp=str2num(char(R.sigasp));
par.viasp=str2num(char(R.viasp));
par.thiasp=str2num(char(R.thiasp));
par.uiasp=str2num(char(R.uiasp));
par.k=str2num(char(R.permea));
par.phi=str2num(char(R.porosity));
par.cf=str2num(char(R.compress));
par.etaf=str2num(char(R.viscosity));
par.rhof=str2num(char(R.rho_f));
par.paraminj=str2num(char(R.paraminj));
par.nx=str2num(char(R.nx));
par.ny=str2num(char(R.ny));
par.dx=str2num(char(R.dx));
par.dy=str2num(char(R.dy));
par.niter=str2num(char(R.niter));
par.nitscreen=str2num(char(R.nit_screen));
par.methinit=str2num(char(R.meth_init));
par.pramdt=str2num(char(R.paramdt));
par.pathinit=R.pathinit;
par.qcat=str2num(char(R.qcat));
par.qmoy=str2num(char(R.qmoy));
par.qprof=str2num(char(R.qprof));
par.qprofhvc=str2num(char(R.qprofhvc));
par.qloc=str2num(char(R.qloc));
par.vfrec=str2num(char(R.vfrec));
par.tfrec=str2num(char(R.tfrec));
par.dtfrec=str2num(char(R.dtfrec));
par.nitrec=str2num(char(R.nitrec));
par.pathres=R.pathres;
par.nploc=str2num(char(R.nploc));
par.xrloc=str2num(char(R.xrloc));
par.yrloc=str2num(char(R.yrloc));
par.hboxm=str2num(char(R.hboxm));
par.pathinjdata=R.pathinjdata;


par.toff=par.paraminj(1);
par.xinj=par.paraminj(2);
par.yinj=par.paraminj(3);
par.qinj=par.paraminj(4);

finfo=ncinfo([rep,'init.nc']);
szcat=size(finfo.Variables);

va=ncread([rep,'init.nc'],'a');
vb=ncread([rep,'init.nc'],'b');
vdc=ncread([rep,'init.nc'],'dc');
vs=ncread([rep,'init.nc'],'s');
vvi=ncread([rep,'init.nc'],'vi');
vthi=ncread([rep,'init.nc'],'thi');
%if szcat(2)==8
%    vpi=ncread([rep,'init.nc'],'pi');
%else
vpi=zeros(size(va));
%end

par.amap=reshape(va,par.ny,par.nx);
par.bmap=reshape(vb,par.ny,par.nx);
par.dcmap=reshape(vdc,par.ny,par.nx);
par.smap=reshape(vs,par.ny,par.nx);
par.vimap=reshape(vvi,par.ny,par.nx);
par.thimap=reshape(vthi,par.ny,par.nx);
par.pimap=reshape(vpi,par.ny,par.nx);

end


%--------------------------------------%
function [t,dt,v,vmax,vmvw,th,tau,u]=extractimeseriesnc(rep)

finfo=ncinfo([rep,'qmoy.nc']);
szcat=size(finfo.Variables);

t=ncread([rep,'qmoy.nc'],'time');
dt=ncread([rep,'qmoy.nc'],'time delai');
v=ncread([rep,'qmoy.nc'],'mean slip rate');
vmax=ncread([rep,'qmoy.nc'],'max slip rate');
th=ncread([rep,'qmoy.nc'],'mean state');
tau=ncread([rep,'qmoy.nc'],'mean shear stress');
u=ncread([rep,'qmoy.nc'],'mean slip');

n=length(t);

if szcat(2)>7
    vmvw=ncread([rep,'qmoy.nc'],'mean slip rate vw');
else
    vmvw=zeros(n,1);
end

end

%--------------------------------------%
function [tloc,vloc,thloc,tauloc,uloc,ploc,qdarcyloc]=extract_local_timeseries_nc(rep,par)

tloc=[];dtloc=[];vloc=[];thloc=[];tauloc=[];uloc=[];ploc=[];qdarcyloc=[];

for i=1:par.nploc
    finfo=ncinfo([rep,'qloc',num2str(i),'.nc']);
    szcat=size(finfo.Variables);

    vt=ncread([rep,'qloc',num2str(i),'.nc'],'time');
    vdt=ncread([rep,'qloc',num2str(i),'.nc'],'time delai');
    vu=ncread([rep,'qloc',num2str(i),'.nc'],'slip');
    vth=ncread([rep,'qloc',num2str(i),'.nc'],'state');
    vv=ncread([rep,'qloc',num2str(i),'.nc'],'slip rate');
    vtau=ncread([rep,'qloc',num2str(i),'.nc'],'shear stress');

    n=length(vt);

    if szcat(2)>6
        vp=ncread([rep,'qloc',num2str(i),'.nc'],'pore pressure');
    else
        vp=zeros(n,1);
    end

    if szcat(2)>7
        vq=ncread([rep,'qloc',num2str(i),'.nc'],'darcy vel');
    else
        vq=zeros(n,1);
    end

    tloc=[tloc vt];
    dtloc=[dtloc vdt];
    vloc=[vloc vv];
    thloc=[thloc vth];
    tauloc=[tauloc vtau];
    uloc=[uloc vu];
    ploc=[ploc vp];
    qdarcyloc=[qdarcyloc vq];

    clear vt vdt vu vth vv vtau vp vq
end

end

%--------------------------------------%
function [t0,dt0,tev,xi,yi,xb,yb,s,m0,dtau,taui,tauf,du]=extractcat(rep,par)


finfo=ncinfo([rep,'earthquake_catalogue.nc']);
szcat=size(finfo.Variables);

t0=ncread([rep,'earthquake_catalogue.nc'],'onset time');
dt0=ncread([rep,'earthquake_catalogue.nc'],'onset time delai');
tev=ncread([rep,'earthquake_catalogue.nc'],'event duration');
xi=ncread([rep,'earthquake_catalogue.nc'],'x initiation');
yi=ncread([rep,'earthquake_catalogue.nc'],'y initiation');
xb=ncread([rep,'earthquake_catalogue.nc'],'x barycenter');
yb=ncread([rep,'earthquake_catalogue.nc'],'y barycenter');
if par.dy*par.dx>0
    s=par.dx*par.dy*ncread([rep,'earthquake_catalogue.nc'],'number of elements');
else
    s=par.dx*par.dx*ncread([rep,'earthquake_catalogue.nc'],'number of elements');
end
m0=ncread([rep,'earthquake_catalogue.nc'],'coseismic moment');
dtau=ncread([rep,'earthquake_catalogue.nc'],'coseismic stress drop');
du=ncread([rep,'earthquake_catalogue.nc'],'coseismic slip');

neqk=length(t0);

if szcat(2)>11
    taui=ncread([rep,'earthquake_catalogue.nc'],'shear stress init');
    tauf=ncread([rep,'earthquake_catalogue.nc'],'shear stress final');
else
    taui=zeros(neqk,1);
    tauf=zeros(neqk,1);
end

end

%-----------------------------------------%

function Result = ini2struct(FileName)
%==========================================================================
%  Author: Andriy Nych ( nych.andriy@gmail.com )
% Version:        733341.4155741782200
%==========================================================================
% 
% INI = ini2struct(FileName)
% 
% This function parses INI file FileName and returns it as a structure with
% section names and keys as fields.
% 
% Sections from INI file are returned as fields of INI structure.
% Each fiels (section of INI file) in turn is structure.
% It's fields are variables from the corresponding section of the INI file.
% 
% If INI file contains "oprhan" variables at the beginning, they will be
% added as fields to INI structure.
% 
% Lines starting with ';' and '#' are ignored (comments).
% 
% See example below for more information.
% 
% Usually, INI files allow to put spaces and numbers in section names
% without restrictions as long as section name is between '[' and ']'.
% It makes people crazy to convert them to valid Matlab variables.
% For this purpose Matlab provides GENVARNAME function, which does
%  "Construct a valid MATLAB variable name from a given candidate".
% See 'help genvarname' for more information.
% 
% The INI2STRUCT function uses the GENVARNAME to convert strange INI
% file string into valid Matlab field names.
% 
% [ test.ini ]~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
%     SectionlessVar1=Oops
%     SectionlessVar2=I did it again ;o)
%     [Application]
%     Title = Cool program
%     LastDir = c:\Far\Far\Away
%     NumberOFSections = 2
%     [1st section]
%     param1 = val1
%     Param 2 = Val 2
%     [Section #2]
%     param1 = val1
%     Param 2 = Val 2
% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% The function converts this INI file it to the following structure:
% 
% [ MatLab session (R2006b) ]~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  >> INI = ini2struct('test.ini');
%  >> disp(INI)
%         sectionlessvar1: 'Oops'
%         sectionlessvar2: 'I did it again ;o)'
%             application: [1x1 struct]
%             x1stSection: [1x1 struct]
%            section0x232: [1x1 struct]
% 
%  >> disp(INI.application)
%                    title: 'Cool program'
%                  lastdir: 'c:\Far\Far\Away'
%         numberofsections: '2'
% 
%  >> disp(INI.x1stSection)
%         param1: 'val1'
%         param2: 'Val 2'
% 
%  >> disp(INI.section0x232)
%         param1: 'val1'
%         param2: 'Val 2'
% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% NOTE.
% WhatToDoWithMyVeryCoolSectionAndVariableNamesInIniFileMyVeryCoolProgramWrites?
% GENVARNAME also does the following:
%   "Any string that exceeds NAMELENGTHMAX is truncated". (doc genvarname)
% Period.
% 
% =========================================================================
Result = [];                            % we have to return something
CurrMainField = '';                     % it will be used later
f = fopen(FileName,'r');                % open file
while ~feof(f)                          % and read until it ends
    s = strtrim(fgetl(f));              % Remove any leading/trailing spaces
    if isempty(s)
        continue;
    end;
    if (s(1)==';')                      % ';' start comment lines
        continue;
    end;
    if (s(1)=='#')                      % '#' start comment lines
        continue;
    end;
    if ( s(1)=='[' ) && (s(end)==']' )
        % We found section
        CurrMainField = genvarname(lower(s(2:end-1)));
        Result.(CurrMainField) = [];    % Create field in Result
    else
        % ??? This is not a section start
        [par,val] = strtok(s, '=');
        val = CleanValue(val);
        if ~isempty(CurrMainField)
            % But we found section before and have to fill it
            Result.(CurrMainField).(lower(genvarname(par))) = val;
        else
            % No sections found before. Orphan value
            Result.(lower(genvarname(par))) = val;
        end
    end
end
fclose(f);
return;

end


function res = CleanValue(s)
res = strtrim(s);
if strcmpi(res(1),'=')
    res(1)=[];
end
res = strtrim(res);
return;

end




