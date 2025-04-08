load('TGdata.mat')



load('TG40.mat')
load('TG50.mat')
load('TG57.mat')

TG40.id=repmat(40,height(TG40),1);
TG50.id=repmat(50,height(TG50),1);
TG57.id=repmat(57,height(TG57),1);

TG40.lat=repmat(TGinfo.Lat(40),height(TG40),1);
TG50.lat=repmat(TGinfo.Lat(50),height(TG50),1);
TG57.lat=repmat(TGinfo.Lat(57),height(TG57),1);


TG40.lon=repmat(TGinfo.Lon(40),height(TG40),1);
TG50.lon=repmat(TGinfo.Lon(50),height(TG50),1);
TG57.lon=repmat(TGinfo.Lon(57),height(TG57),1);


tg=[TG40;TG50;TG57];
clearvars -except tg  
tg.t=decyear(tg.Time);
tg.Properties.VariableNames{3} = 'lu';



load('F:\Analysis\BS\LU\ers2\588.mat')
said(1:height(TBL),1)=1;
TBL = addvars(TBL,said);
sa=TBL;
clear said
load('F:\Analysis\BS\LU\envisat\588.mat')
said(1:height(TBL),1)=2;
TBL = addvars(TBL,said);
sa=[sa;TBL];
clear said
load('F:\Analysis\BS\LU\saral\588.mat')
said(1:height(TBL),1)=3;
TBL = addvars(TBL,said);
sa=[sa;TBL];
clear said
load('F:\Analysis\BS\LU\s3a\158.mat')
said(1:height(TBL),1)=4;
TBL = addvars(TBL,said);
sa=[sa;TBL];
clear said
clear h

sa.date=datetime(datevec(ConvertSerialYearToDate(sa.t)));
sa.y=year(sa.date);

sa = rmmissing(sa,'DataVariables',{'ssh'}); %remove NaN data

clear TBL
%% 
DEG=table();
for j=1:4
id=sa.id;
id=unique(id);
for k=id(1):id(end)
       H = height(DEG);
       
       if~isempty(unique(sa.pass(sa.id==k&sa.said==j)))
           DEG.id(H+1)=j;
           DEG.cycle(H+1)=unique(sa.cycle(sa.id==k&sa.said==j));
           time=sa.date(sa.id==k&sa.said==j);
           time=time(1);
           DEG.month(H+1)=month(time);
           DEG.year(H+1)=year(time);
           
           track_selected=[sa.lat(sa.id==k&sa.said==j),...
                           (sa.ssh(sa.id==k&sa.said==j)-sa.N(sa.id==k&sa.said==j)+sa.Ell_corr(sa.id==k&sa.said==j)+sa.dac(sa.id==k&sa.said==j))*100,...
                           sa.t(sa.id==k&sa.said==j)];
        
        if ~isempty(track_selected(~any(ismissing(track_selected(:,2)),2),:))
           
            
            %   remove gross errors
            [m1,~]=find(abs(track_selected(:,2))>=200);
            track_selected(m1,:)=NaN;
            
            %   remove random errors
            [m11,~]=find(abs(track_selected(:,2))>=mean(abs(track_selected(:,2)),'omitnan')+std(abs(track_selected(:,2)),'omitnan')*3);
            track_selected(m11,:)=NaN;
            
            %   remove outlier data using scaled MAD
            [~,out]=rmoutliers(track_selected(:,2),'movmedian',100);
            [m2, ~]=find(out==1);
            track_selected(m2,:)=NaN;
            
            [m3,~]=find(isnan(track_selected(:,1)));
            track_selected(m3,:)=[];
            
            
            tg_lu=griddata(tg.lat,tg.t,tg.lu,track_selected(:,1),track_selected(:,3));
            tg_obs=griddata(tg.lat,tg.t,tg.obs,track_selected(:,1),track_selected(:,3));

            
            [m,~]=find(isnan(tg_lu));
            tg_obs(m,:)=[];
            tg_lu(m,:)=[];
            track_selected(m,:)=[];
            
            dh_lu=track_selected(:,2)-tg_lu;
            p1 = polyfit(track_selected(:,1),dh_lu,1);
            f1 = polyval(p1,track_selected(:,1));
            
            ff2=fitlm(track_selected(:,1),dh_lu);
            deltadh=predict(ff2,track_selected(end,1))-predict(ff2,track_selected(1,1));
            deltaphi=track_selected(end,1)-track_selected(1,1)*110*1000*100;
            DEG.lu(H+1)=rad2deg(atan(deltadh/deltaphi))*60*60*60*100;

            
            dh_obs=track_selected(:,2)-tg_obs;
            p2 = polyfit(track_selected(:,1),dh_obs,1);
            f2 = polyval(p2,track_selected(:,1));
            ff1=fitlm(track_selected(:,1),dh_obs);
            deltadh=predict(ff1,track_selected(end,1))-predict(ff1,track_selected(1,1));
            deltaphi=track_selected(end,1)-track_selected(1,1)*110*1000*100;
            DEG.obs(H+1)=rad2deg(atan(deltadh/deltaphi))*60*60*60*100;
%             
%             plot(track_selected(:,1),dh_obs,'.','Color',[.2 1 .2])
%             hold on
%             plot(track_selected(:,1),f2,'-','Color',[0.4660 0.6740 0.1880],'LineWidth',2.5)
%             
%             plot(track_selected(:,1),dh_lu,'.','Color',[1 .5 .5])
%             plot(track_selected(:,1),f1,'-r','LineWidth',2.5)
%            
%             
%             
            

clear ff1
clear ff2
clear f1
clear f2
clear p1
clear p2
clear dh_lu
clear dh_obs
clear deltaphi
clear deltadh
clear m
clear m1
clear m11
clear m2
clear m3
clear track_selected
clear tg_lu
clear tg_obs
clear out

        end
       end
end
end
%% plot


plot(decyear(DEG.time),DEG.lu/10,'.r')
hold on
plot(decyear(DEG.time),DEG.obs/10,'.b')

p1 = polyfit(decyear(DEG.time),DEG.lu/10,1);
f1 = polyval(p1,decyear(DEG.time));

p2 = polyfit(decyear(DEG.time),DEG.obs/10,1);
f2 = polyval(p2,decyear(DEG.time));

plot(decyear(DEG.time),f1,'-r','LineWidth',2.5)
plot(decyear(DEG.time),f2,'-b','LineWidth',2.5)
yline(0,'--','Refrence','LineWidth',1.5,'Color',[.3 .3 .3])
xline(2000,'-','epoch2000','LineWidth',1.5,'Color',[.3 .3 .3])     
xticks(1995:3:2020)
set(gca,'XGrid','on')

text(2013,-10,strcat(sprintf('y = %0.2f x + %0.2f',p2(1),p2(2))),'Color','b','FontName','Times New Roman','FontWeight','Bold','FontSize',13)
text(2013,-15,strcat(sprintf('y = %0.2f x + %0.2f',p1(1),p1(2))),'Color','r','FontName','Times New Roman','FontWeight','Bold','FontSize',13)



ylabel('Tilt [" \times 10]','FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;



 