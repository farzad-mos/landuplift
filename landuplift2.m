% load TG data

load('TG40.mat')
load('TG50.mat')
load('TG51.mat')
load('TG552.mat')
load('TG55.mat')
load('TG45.mat')
load('TG57.mat')


load('TGdata.mat')

%% plot SA pass
load('F:\Analysis\BS\LU\saral\588.mat')

k=1;
geoplot(TBL.lat(TBL.id==k),TBL.lon(TBL.id==k),'.b')
hold on
load('F:\Analysis\BS\LU\s3a\158.mat')
geoplot(TBL.lat(TBL.id==k),TBL.lon(TBL.id==k),'.', 'Color',[0.3010 0.7450 0.9330])
ID=[40,50,57];
for i=1:3
geoplot(TGinfo.Lat(ID(i)),TGinfo.Lon(ID(i)),'^k','MarkerFaceColor','k','MarkerSize',10)
text(TGinfo.Lat(ID(i)),TGinfo.Lon(ID(i))-1.2,num2str(ID(i)),'Color','r','FontSize',14)
end

ax=gca; ax.FontSize=14; ax.FontWeight='Bold'; grid on
set(gca,'fontname','Times New Roman');



%% plot TG data
t1=datetime(1995,1,1,1,0,0);
t2=datetime(2020,1,1,1,0,0);

subplot(3,1,1)
plot(TG40.Time,TG40.obs,'DisplayName','TG40.obs');hold on;plot(TG40.Time,TG40.LU,'DisplayName','TG40.LU')
xlim([t1 t2])
% xticks(t1:3:t2)
set(gca,'XGrid','on')
ylabel('DT [cm]','FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;

subplot(3,1,2)
plot(TG50.Time,TG50.obs,'DisplayName','TG50.obs');hold on;plot(TG50.Time,TG50.LU,'DisplayName','TG50.LU')
xlim([t1 t2])
% xticks(t1:3:t2)
set(gca,'XGrid','on')
ylabel('DT [cm]','FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;

subplot(3,1,3)
plot(TG57.Time,TG57.obs,'DisplayName','TG57.obs');hold on;plot(TG57.Time,TG57.LU,'DisplayName','TG57.LU')
xlim([t1 t2])
ylabel('DT [cm]','FontSize',18,'FontWeight','bold');

% xticks(t1:3:t2)
set(gca,'XGrid','on')
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;


% 
% subplot(7,1,4)
% plot(TG51.Time,TG51.obs,'DisplayName','TG51.obs');hold on;plot(TG51.Time,TG51.LU,'DisplayName','TG51.LU')
% xlim([t1 t2])
% 
% subplot(7,1,5)
% plot(TG55.Time,TG55.obs,'DisplayName','TG55.obs');hold on;plot(TG55.Time,TG55.LU,'DisplayName','TG55.LU')
% xlim([t1 t2])
% 
% subplot(7,1,6)
% plot(TG552.Time,TG552.obs,'DisplayName','TG552.obs');hold on;plot(TG552.Time,TG552.LU,'DisplayName','TG552.LU')
% xlim([t1 t2])
% 
% subplot(7,1,7)
% plot(TG57.Time,TG57.obs,'DisplayName','TG57.obs');hold on;plot(TG57.Time,TG57.LU,'DisplayName','TG57.LU')
% xlim([t1 t2])

%% select TG ID
TGID=55;

if TGID==552
%Ratan TG only
Lon=20.8950;
Lat=63.9861;
else
Lon=TGinfo.Lon(TGID);
Lat=TGinfo.Lat(TGID);
end

TG55.y=year(TG55.Time);

%% extract nearby data

load('F:\Analysis\BS\LU\ers2\588.mat')
[m,~]=find(TBL.lat<64);
TBL(m,:)=[];

for i=1:height(TBL)
TBL.dist(i)=distance([TBL.lat(i) TBL.lon(i)],[Lat Lon],referenceEllipsoid('WGS84'))/1000; %measure TG to TBL points diTBLnce
end
TBL=TBL(TBL.dist<=70,:);
ers=TBL;


load('F:\Analysis\BS\LU\envisat\588.mat')
[m,~]=find(TBL.lat<64);
TBL(m,:)=[];
for i=1:height(TBL)
TBL.dist(i)=distance([TBL.lat(i) TBL.lon(i)],[Lat Lon],referenceEllipsoid('WGS84'))/1000; %measure TG to TBL points diTBLnce
end

TBL=TBL(TBL.dist<=70,:);
envi=TBL;


load('F:\Analysis\BS\LU\saral\588.mat')
[m,~]=find(TBL.lat<64);
TBL(m,:)=[];
for i=1:height(TBL)
TBL.dist(i)=distance([TBL.lat(i) TBL.lon(i)],[Lat Lon],referenceEllipsoid('WGS84'))/1000; %measure TG to TBL points diTBLnce
end
TBL=TBL(TBL.dist<=70,:);
saral=TBL;


load('F:\Analysis\BS\LU\s3a\158.mat')
[m,~]=find(TBL.lat<64);
TBL(m,:)=[];
for i=1:height(TBL)
TBL.dist(i)=distance([TBL.lat(i) TBL.lon(i)],[Lat Lon],referenceEllipsoid('WGS84'))/1000; %measure TG to TBL points diTBLnce
end
TBL=TBL(TBL.dist<=70,:);
s3a=TBL;

clear TBL
%% 
sa=table();
for i=1:4
    
    if i==1
        SA_cor=ers;
    elseif i==2
        SA_cor=envi;
        
    elseif i==3
        SA_cor=saral;
        
    elseif i==4
        SA_cor=s3a;
        
    end
    
    
    id=SA_cor.id;
    id=unique(id);
    
    for k=id(1):id(end)
        H = height(sa);
        
        sa.id(H+1)=i;
        sa.cycle(H+1)=unique(SA_cor.cycle(SA_cor.id==k));
        time=SA_cor.t(SA_cor.id==k);
        sa.t(H+1)=time(1);
        
        track_selected=[(SA_cor.ssh(SA_cor.id==k)-SA_cor.N(SA_cor.id==k)+SA_cor.Ell_corr(SA_cor.id==k)+SA_cor.dac(SA_cor.id==k))*100,...
                        (SA_cor.ssh(SA_cor.id==k)-SA_cor.N(SA_cor.id==k)+SA_cor.Ell_corr(SA_cor.id==k))*100];
        %   remove gross errors
        [m1,~]=find(abs(track_selected(:,1))>=100);
        track_selected(m1,:)=NaN;
        
        %   remove random errors
        [m11,~]=find(abs(track_selected(:,1))>=mean(abs(track_selected(:,1)),'omitnan')+std(abs(track_selected(:,1)),'omitnan')*3);
        track_selected(m11,:)=NaN;
        
        %   remove outlier data using scaled MAD
        [~,out]=rmoutliers(track_selected(:,1),'movmedian',200);
        [m2, ~]=find(out==1);
        track_selected(m2,:)=NaN;
     
     
     sa.dt(H+1)=mean(track_selected(:,1),'omitnan');
     sa.dtdac(H+1)=mean(track_selected(:,2),'omitnan');
     
     clear time
     clear track_selected
     clear m1
     clear m11
     clear m2
     
    end
     
end
     
sa.date= datetime(datevec(ConvertSerialYearToDate(sa.t)));
sa.tglu=interp1(decyear(TG57.Time),TG57.LU,sa.t);
sa.tgobs=interp1(decyear(TG57.Time),TG57.obs,sa.t);

sa = rmmissing(sa,'DataVariables',{'dt'}); %remove NaN data
sa.y=year(sa.date);

%% dt
p1 = polyfit(sa.t,sa.dt,1);
f1 = polyval(p1,sa.t);
            
p2 = polyfit(sa.t,sa.tglu,1);
f2 = polyval(p2,sa.t);

p3 = polyfit(sa.t,sa.tgobs,1);
f3 = polyval(p3,sa.t);
            
plot(sa.t,sa.dt,'.k');
hold on
plot(sa.t,sa.tglu,'.r');
plot(sa.t,sa.tgobs,'.b')

plot(sa.t,f1,'-k','LineWidth',2.5)
plot(sa.t,f2,'-r','LineWidth',2.5)
plot(sa.t,f3,'-b','LineWidth',2.5)





%% diff

for i=1:3
    if i==1
        load('satg40.mat')
    elseif i==2
        load('satg50.mat')
    elseif i==3
        load('satg57.mat')
    end
    
    subplot(3,1,i)
    dh_lu=sa.dt-sa.tglu;
    dh_obs=sa.dt-sa.tgobs;
    
    p1 = polyfit(sa.t,dh_lu,1);
    f1 = polyval(p1,sa.t);
    
    p2 = polyfit(sa.t,dh_obs,1);
    f2 = polyval(p2,sa.t);
    plot(sa.t,dh_lu,'.r');
    hold on
    plot(sa.t,dh_obs,'.b');
    
    plot(sa.t,f1,'-r','LineWidth',2.5)
    plot(sa.t,f2,'-b','LineWidth',2.5)
    yline(0,'-','Refrence','LineWidth',1.5,'Color',[.3 .3 .3])
    xline(2000,'-','epoch','LineWidth',1.5,'Color',[.3 .3 .3])
    ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
    ylabel('\Delta DT [cm]','FontSize',18,'FontWeight','bold');
    if i~=3
        set(gca,'XGrid','on','xticklabel',[])
    else
        set(gca,'XGrid','on')
    end
    
    xticks(1995:2020)
    ylim([-80 80])
    hold off
    
    clear sa
    
end

%% msl
clear TG
clear msl

sumer=0;

% TG=TG40;
% load('satg40.mat')
% TG=TG50;
% load('satg50.mat')
TG=TG57;
load('satg57.mat')

if sumer==1
% REMOVE SUMMER DAYS
TG.m=month(TG.Time);
sa.m=month(sa.date);
sa.dt(~ismember(sa.m,[6,7,8]))=nan;
TG.LU(~ismember(TG.m,[6,7,8]))=nan;
sa = rmmissing(sa,'DataVariables',{'dt'}); %remove NaN data
TG = rmmissing(TG,'DataVariables',{'LU'}); %remove NaN data
end



[group,year]=findgroups(sa.y);
samsl=splitapply(@mean, sa.dt, group);
samsldac=splitapply(@mean, sa.dtdac, group);


[group,year2]=findgroups(TG.y);
tgmsl=[splitapply(@mean, TG.LU, group),splitapply(@mean, TG.obs, group)];

lumsl=interp1(year2,tgmsl(:,1),year);
obsmsl=interp1(year2,tgmsl(:,2),year);


msl=table(year,samsl,samsldac,lumsl,obsmsl);


%% plot MSL



% decorrected DAC
figure(1)
subplot(2,1,1)     
plot(msl.year,msl.samsl,'-.k');
hold on
plot(msl.year,msl.lumsl,'-.r')
plot(msl.year,msl.obsmsl,'-.b')
f1=fitlm(msl.year,msl.lumsl);
f2=fitlm(msl.year,msl.obsmsl);
f3=fitlm(msl.year,msl.samsl);
h1=plot(f1);
h2=plot(f2);
h3=plot(f3);

fitHandle = findobj(h1,'DisplayName','Fit');
fitHandle.Color ='r';
h1(end-1,1).Visible='off';
h1(end,1).Visible='off';
dataHandle = findobj(h1,'DisplayName','Data');
dataHandle.Color = 'r'; 

fitHandle = findobj(h2,'DisplayName','Fit');
fitHandle.Color ='b';
h2(end-1,1).Visible='off';
h2(end,1).Visible='off';
dataHandle = findobj(h2,'DisplayName','Data');
dataHandle.Color = 'b'; 

fitHandle = findobj(h3,'DisplayName','Fit');
fitHandle.Color ='k';
h3(end-1,1).Visible='off';
h3(end,1).Visible='off';
dataHandle = findobj(h3,'DisplayName','Data');
dataHandle.Color = 'k'; 

p(1)=plot(msl.year,msl.samsl,'.k','DisplayName',strcat('SA trend= ',num2str(f3.Coefficients.Estimate(2)*10,2),' mm/year'));
p(2)=plot(msl.year,msl.lumsl,'.r','DisplayName',strcat('TG LU trend= ',num2str(f1.Coefficients.Estimate(2)*10,2),' mm/year'));
p(3)=plot(msl.year,msl.obsmsl,'.b','DisplayName',strcat('TG obs trend= ',num2str(f2.Coefficients.Estimate(2)*10,2),' mm/year'));

xticks(1995:3:2020)
set(gca,'XGrid','on')
xline(2000,'-','epoch2000','LineWidth',1.5,'Color',[.3 .3 .3])     
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
ylabel('MDT [cm]','FontSize',18,'FontWeight','bold');
ylim([-10 40])
legend(p(1:3),'FontSize',10)
xlim([1993 2021])

subplot(2,1,2)     
plot(msl.year,msl.samsl-msl.lumsl,'-.r','LineWidth',1.5)
hold on
plot(msl.year,msl.samsl-msl.obsmsl,'-.b','LineWidth',1.5)

yline(0,'-','Refrence','LineWidth',1.5,'Color',[.3 .3 .3])
xline(2000,'-','epoch2000','LineWidth',1.5,'Color',[.3 .3 .3])     
xticks(1995:3:2020)
set(gca,'XGrid','on')
ylim([-20 20])
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;     
ylabel('\Delta MDT [cm]','FontSize',18,'FontWeight','bold');
xlim([1993 2021])     
rmse1=rms(msl.samsldac-msl.lumsl);
rmse2=rms(msl.samsldac-msl.obsmsl);
text(1994,15,strcat('RMSE= ',num2str(rmse1,3),' cm'),'Color','r','FontSize',14,'FontWeight','bold')
text(1994,10,strcat('RMSE= ',num2str(rmse2,3),' cm'),'Color','b','FontSize',14,'FontWeight','bold')
     
bias_lu=predict(f1,2000)-predict(f3,2000)     
bias_obs=predict(f2,2000)-predict(f3,2000)     


% with DAC
figure(2)
subplot(2,1,1)     
plot(msl.year,msl.samsldac,'-.k');
hold on
plot(msl.year,msl.lumsl,'-.r')
plot(msl.year,msl.obsmsl,'-.b')

f1=fitlm(msl.year,msl.lumsl);
f2=fitlm(msl.year,msl.obsmsl);
f3=fitlm(msl.year,msl.samsldac);


h1=plot(f1);
h2=plot(f2);
h3=plot(f3);


fitHandle = findobj(h1,'DisplayName','Fit');
fitHandle.Color ='r';
h1(end-1,1).Visible='off';
h1(end,1).Visible='off';
dataHandle = findobj(h1,'DisplayName','Data');
dataHandle.Color = 'r'; 

fitHandle = findobj(h2,'DisplayName','Fit');
fitHandle.Color ='b';
h2(end-1,1).Visible='off';
h2(end,1).Visible='off';
dataHandle = findobj(h2,'DisplayName','Data');
dataHandle.Color = 'b'; 

fitHandle = findobj(h3,'DisplayName','Fit');
fitHandle.Color ='k';
h3(end-1,1).Visible='off';
h3(end,1).Visible='off';
dataHandle = findobj(h3,'DisplayName','Data');
dataHandle.Color = 'k'; 

p(1)=plot(msl.year,msl.samsldac,'.k','DisplayName',strcat('SA_D_A_C trend= ',num2str(f3.Coefficients.Estimate(2)*10,2),' mm/year'));
p(2)=plot(msl.year,msl.lumsl,'.r','DisplayName',strcat('TG LU trend= ',num2str(f1.Coefficients.Estimate(2)*10,2),' mm/year'));
p(3)=plot(msl.year,msl.obsmsl,'.b','DisplayName',strcat('TG obs trend= ',num2str(f2.Coefficients.Estimate(2)*10,2),' mm/year'));

xticks(1995:3:2020)
set(gca,'XGrid','on')
xline(2000,'-','epoch2000','LineWidth',1.5,'Color',[.3 .3 .3])     
ylim([-10 40])
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
ylabel('MDT [cm]','FontSize',18,'FontWeight','bold');
xlim([1993 2021])
legend(p(1:3),'FontSize',10)



subplot(2,1,2)     
plot(msl.year,msl.samsldac-msl.lumsl,'-.r','LineWidth',1.5)
hold on
plot(msl.year,msl.samsldac-msl.obsmsl,'-.b','LineWidth',1.5)

p1 = polyfit(msl.year,msl.samsldac-msl.lumsl,1);
ff1 = polyval(p1,msl.year);

p2 = polyfit(msl.year,msl.samsldac-msl.obsmsl,1);
ff2 = polyval(p2,msl.year);

plot(msl.year,ff1,'-r','LineWidth',2.5)
plot(msl.year,ff2,'-b','LineWidth',2.5)

yline(0,'-','Refrence','LineWidth',1.5,'Color',[.3 .3 .3])
xline(2000,'-','epoch2000','LineWidth',1.5,'Color',[.3 .3 .3])     
xticks(1995:3:2020)
set(gca,'XGrid','on')
xlim([1993 2021])
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;     
ylabel('\Delta MDT [cm]','FontSize',18,'FontWeight','bold');
ylim([-30 30])
rmse1=rms(msl.samsldac-msl.lumsl);
rmse2=rms(msl.samsldac-msl.obsmsl);
text(1994,15,strcat('RMSE= ',num2str(rmse1,3),' cm'),'Color','r','FontSize',14,'FontWeight','bold')
text(1994,10,strcat('RMSE= ',num2str(rmse2,3),' cm'),'Color','b','FontSize',14,'FontWeight','bold')

bias_lu_dac=predict(f1,2000)-predict(f3,2000)     
bias_obs_dac=predict(f2,2000)-predict(f3,2000)     
