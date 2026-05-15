clear;
%close all;
dim=2;
dimp1=dim+1;
rx3=1/sqrt(5);
% sphere reflecting surface
h=1;
hr=2;
f=1;
g=@(z)f*h./(h^2+z(:,1).^2+z(:,2).^2).^1.5;
rtar=h*sqrt((1-rx3^2)/rx3^2);
%%
left=[-rtar -rtar];
right=[rtar rtar];
Nm=17;%3 5 8 10 13 16 20 25 28
coord1=linspace(left(1),right(1),Nm);
coord2=linspace(left(2),right(2),Nm);
[X,Y]=meshgrid(coord1,coord2);
Z=[X(:) Y(:)];
Z=Z(sqrt(sum(Z.^2,2))<=rtar,:);
den=g(Z);

% nZ=3;%13 40 73 142 218 349 559 561 712
% H = Hammersley(nZ, dim)';
% %figure
% %plot(H(1,:),H(2,:),'ro')
% %H=sobol(1:K1,:)';
% Z=-[rtar rtar]+2*rtar*H;
% Z=Z(sqrt(sum(Z.^2,2))<=rtar,:);
% den=g(Z);

% nZ=3;%16 41 75 138(250) 220 355 600
% Z=-[rtar rtar]+2*rtar*rand(nZ,dim);
% %[~,Z] = kmeans(Z,nZ/8);
% Z=Z(sqrt(sum(Z.^2,2))<=rtar,:);
% den=g(Z);

% Z=[0 -1;0 1];
% den=g(Z);

den1=den/sum(den);
errorbound=1e-4;
ind=den1>errorbound;
target_den=den(ind);
K=length(target_den);
x3=-h;
F2=[Z(ind,:) x3*ones(K,1)];
%plot3(F2(:,1),F2(:,2),target_den,'ro')
figure
plot(F2(:,1),F2(:,2),'ro')
%%
[target_den,sindtarget_den]=sort(target_den,'descend');%'ascend''descend'
F2=F2(sindtarget_den(1:K),:);

Sample_num=23e4;%60
dray=sample_on_hemisphere(Sample_num,dimp1);
dray=dray(:,dray(3,:)>=rx3);

F1=[0 0 0];
normF2=sqrt(sum(F2.^2,2))';
normalizeF2=F2./normF2';
M=max(normF2);

gami=max(dray'*normalizeF2');
gam=max(gami);

a1=sqrt((2*hr+h)^2+F2(1,1)^2+F2(1,2)^2)/2;
c1=norm(F2(1,:))/2;
d1=a1-c1^2/a1;

amax=2/(1-gami(1)*(sqrt(1+d1^2/M^2)-d1/M));
amin=(1-gam)/2;
%sphere reflector
d=[d1;amax*d1*ones(K-1,1)];
de=[0;ones(K-1,1)]*(d1*(amax-amin)/4);

I=1*ones(length(dray),1);
I=I*sum(target_den)/sum(I);
errnorm=sum(I);
eplison=errorbound;
eplison2=errorbound*1e-1;

v=F2-F1;
dv=sqrt(sum(v.^2,2));
nv=v./dv;

d=d-de;
C=20e3;
interval_target=target_den+eplison/sqrt(K*(K-1));
tic
for k=1:C
    %modified Kochengrin algorithm
    e=sqrt(1+d.^2./dv.^2)-d./dv;
    t=d./(1-e.*nv*dray);
    mint_ind=t==min(t);
    color=sum(mint_ind.*I',2);
%    color=sum(t==min(t),2)*I; %sphere reflector
%     color=sum(t==max(t),2)*bright*I;%falt reflector
    err=sqrt(sum((color-target_den).^2)/K)/errnorm;
    if err<eplison || norm(de)<eplison2
        break
    end
    lo=color>interval_target;
    lo(1)=0;
    de(lo)=de(lo)/2;
    d(lo)=d(lo)+de(lo);
    if max(lo)
        continue
    end
    lo=color<interval_target;
    lo(1)=0;
    de(lo)=1.25*de(lo);
    d(lo)=d(lo)-de(lo);
end
RefTime=toc

%% sampling
% num_sample_refl=20e4;
% refrays=sample_on_hemisphere(num_sample_refl,dimp1);
% refrays=refrays(:,refrays(3,:)>=rx3)';
% refrays=refrays(1:10e4,:);

num_sample_refl=5e4;%length(refrays);

F1=[0 0 0];
v=F2-F1;
dv=sqrt(sum(v.^2,2));
nv=v./dv;
e=sqrt(1+d.^2./dv.^2)-d./dv;

eps=1e-6;%use:1e-6;
rho=@(x)-eps*lse(-(d./(1-sum(e.*nv.*x,2)))/eps);
gradrho=@(x)gradfg(x,nv,e,d,dimp1,eps);
refdri=zeros(num_sample_refl,dimp1);
HitTar=zeros(num_sample_refl,dimp1);
rhox=zeros(num_sample_refl,1);
for i=1:num_sample_refl
    Drho=gradrho(refrays(i,:));
    Z(i,:)=Drho;
    rhox(i)=rho(refrays(i,:));
    normdir=[Drho(1:dim) 0]-refrays(i,:)*(rhox(i)+Drho(1:dim)*refrays(i,1:dim)');
    normdir=normdir/sqrt(sum(normdir.^2));
    refdri(i,:)=refrays(i,:)-2*(refrays(i,:)*normdir')* normdir;
    distan=(x3-refrays(i,dimp1).*rhox(i))./refdri(i,dimp1);
    HitTar(i,:)=refrays(i,:)*rhox(i)+refdri(i,:)*distan;
end
figure
plot3(refdri(:,1),refdri(:,2),refdri(:,3),'r.')
figure
plot(HitTar(:,1),HitTar(:,2),'r.')

% HitTar1=-h*refrays./refrays(:,dimp1);
% figure
% plot(HitTar1(:,1),HitTar1(:,2),'r.')

%%
Nps=50;
U=linspace(-pi/2,pi/2,Nps);
V=linspace(0,2*pi,Nps);
[sU,sV]=meshgrid(U,V);
sX=2*cos(sU).*cos(sV);
sY=2*cos(sU).*sin(sV);
sZ=2*sin(sU);
figure
surf(sX,sY,sZ,'FaceColor','red','FaceAlpha',0.1,'EdgeColor','none')
hold on

U=linspace(asin(rx3),pi/2,Nps);
V=linspace(0,2*pi,Nps);
[sU,sV]=meshgrid(U,V);
sX=cos(sU).*cos(sV);
sY=cos(sU).*sin(sV);
sZ=sin(sU);
refrays1=[sX(:) sY(:) sZ(:)];
num_sample_refl=length(refrays1);
rhox1=zeros(num_sample_refl,1);
for i=1:num_sample_refl   
    rhox1(i)=rho(refrays1(i,:));
end
refrays1rho=refrays1.*rhox1;
sXp=reshape(refrays1rho(:,1), size(sX));
sYp=reshape(refrays1rho(:,2), size(sX));
sZp=reshape(refrays1rho(:,3), size(sX));
% figure
surf(sXp,sYp,sZp,'EdgeColor','none')
%colormap(mycolormap('viridis'))
camlight right;
 lighting phong;
% shading interp
axis('equal', 'tight')
set(gca,'box','off');
set(gca,'visible','off');

%%
exportgraphics(gca,['C:\Users\Administrator\OneDrive\文档\W2\ReflectorBayesia_v6\tuii\'.....
    'testexa\testexaRefLowEe-4K561.png'],'ContentType','image','Resolution',150)

%%
function res = lse(v)
    % Compute log sum exp
    maxv= max(v); % Max for numerical stability
    res =maxv + log(sum(exp(v - maxv))); % Log-sum-exp operation
end
function res = gradfg(x,nv,e,d,dimp1,eps)
v=-(d./(1-sum(e.*nv.*x,2)))/eps;
maxv=max(v);
vtmax=v-maxv;
expvtmaxden=exp(vtmax);
int1=sum(expvtmaxden);
int2=sum(expvtmaxden.*(d.*e.*(nv-x.*nv(:,dimp1)/x(dimp1))./(1-e.*sum(x.*nv,2)).^2));
res=int2/int1;
end
function [randd,nd]=sample_on_hemisphere(N,dim)
randd=randn(N,dim);
nd=sqrt(sum(randd.^2,2));
randd=randd./nd;
randd(:,3)=abs(randd(:,3));
randd=randd.';
end