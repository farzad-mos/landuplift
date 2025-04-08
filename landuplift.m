load('Z:\data\processed_data\TG data with Land-Uplift correction\TGdata.mat')
load('ProcessedTGdata.mat') %load pre-processed TG records (simulated data)
TG(:,33)=TG(:,33)+15; %Bornholm tide gauge vertical datum differs from the present (common BSCD) one by 15 cm

load('272.mat')
load('158.mat')
load('386.mat')
load('F:\Analysis\BS\finalsa\tg_s3a.mat')
tg=tg_s3a;

%% 


clear tg_selected
clear tg_lu

TBL = rmmissing(TBL,'DataVariables',{'ssh'}); %remove NaN data
id=TBL.id;
id=unique(id);
TBL.date= datetime(datevec(ConvertSerialYearToDate(TBL.t)));

SA_cor=TBL;
p=SA_cor.pass(1); %SA pass no
[n1,~]=find(tg(:,1)==p); %find tg stations


%% RMSE

RMSE=table(); % RMSE table

for k=id(1):id(end)
     H = height(RMSE);
          
     
    if~isempty(unique(SA_cor.pass(SA_cor.id==k)))
        cycle=unique(SA_cor.cycle(SA_cor.id==k));
        time=SA_cor.date(SA_cor.id==k);
        time=time(1);
        RMSE.month(H+1)=month(time);
        RMSE.year(H+1)=year((time));
        str=datestr(time,'mmm yyyy');
        
        % Lat SA HDM
        track_selected=[SA_cor.lat(SA_cor.id==k),(SA_cor.ssh(SA_cor.id==k)-SA_cor.N(SA_cor.id==k)+SA_cor.Ell_corr(SA_cor.id==k))*100,...
            (SA_cor.ssh(SA_cor.id==k)-SA_cor.N(SA_cor.id==k)+SA_cor.Ell_corr(SA_cor.id==k)+SA_cor.dac(SA_cor.id==k))*100,...
            (SA_cor.ssh(SA_cor.id==k)-SA_cor.N(SA_cor.id==k)+SA_cor.Ell_corr(SA_cor.id==k)-SA_cor.oct(SA_cor.id==k))*100,...
            (SA_cor.ssh(SA_cor.id==k)-SA_cor.N(SA_cor.id==k)+SA_cor.Ell_corr(SA_cor.id==k)+SA_cor.dac(SA_cor.id==k)-SA_cor.oct(SA_cor.id==k))*100];
        
        if ~isempty(track_selected(~any(ismissing(track_selected(:,3)),2),:))
            
            t=SA_cor.t(SA_cor.id==k);
            t=t(1);
            
            y=2;x=1;
            while tg(n1,y)~=0
                tg_lu(x,:)=[TGinfo.Lat(tg(n1,y)),TGinfo.Lon(tg(n1,y)),interp1(decyear(TGdate),TG_wLU(:,tg(n1,y)),t),tg(n1,y)];
                tg_selected(x,:)=[TGinfo.Lat(tg(n1,y)),TGinfo.Lon(tg(n1,y)),interp1(decyear(TGdate),TG(:,tg(n1,y)),t),tg(n1,y)];
                x=x+1;
                y=y+1;
            end
            
            tg_lu = sortrows(tg_lu,1);
            tg_selected = sortrows(tg_selected,1);
            
            
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
            
            
            figure(k)
            subplot(3,1,1)
            plot(track_selected(:,1),track_selected(:,2),'.','Color',[0.3010 0.7450 0.9330])
            hold on
            plot(track_selected(:,1),smoothdata(track_selected(:,2),'gaussian',50),'.b')

                        
            %with LU
            plot(tg_lu(:,1),tg_lu(:,3),'or')
            plot(tg_lu(:,1),tg_lu(:,3),'--r')
            text(tg_lu(:,1),tg_lu(:,3)+5,num2str(tg_lu(:,4)),'Color','r','FontSize',12,'FontWeight','Bold')
            
            %withouth LU
            plot(tg_selected(:,1),tg_selected(:,3),'o','Color',[0.4660 0.6740 0.1880])
            plot(tg_selected(:,1),tg_selected(:,3),'--','Color',[0.4660 0.6740 0.1880])
            ylabel('DT [cm]','FontSize',18,'FontWeight','bold');
            ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
            rmse1=rms(track_selected(:,2)-(interp1(tg_lu(:,1),tg_lu(:,3),track_selected(:,1))),'omitnan');
            rmse2=rms(track_selected(:,2)-(interp1(tg_selected(:,1),tg_selected(:,3),track_selected(:,1))),'omitnan');
            text(tg_lu(1,1),tg_lu(1,3)+50,strcat('RMSE[cm]=',num2str(rmse1)),'Color','r','FontSize',14)
            text(tg_lu(1,1),tg_lu(1,3)+30,strcat('RMSE[cm]=',num2str(rmse2)),'Color',[0.4660 0.6740 0.1880],'FontSize',14)
            yyaxis right
            ylabel('+ DAC','Color','k')
            set(gca,'yticklabel',[],'YColor','k')
            title(strcat('Cycle#',num2str(cycle),': ',str))           
            hold off
            
            RMSE.daclu(H+1)=rmse1;
            RMSE.dacw(H+1)=rmse2;
            clear rmse1
            clear rmse2
            
            
            %with dac
            subplot(3,1,2)
            plot(track_selected(:,1),track_selected(:,3),'.','Color',[0.3010 0.7450 0.9330])
            hold on
            plot(track_selected(:,1),smoothdata(track_selected(:,3),'gaussian',50),'.b')
            %with LU
            plot(tg_lu(:,1),tg_lu(:,3),'or')
            plot(tg_lu(:,1),tg_lu(:,3),'--r')
            %withouth LU
            plot(tg_selected(:,1),tg_selected(:,3),'o','Color',[0.4660 0.6740 0.1880])
            plot(tg_selected(:,1),tg_selected(:,3),'--','Color',[0.4660 0.6740 0.1880])
            ylabel('DT [cm]','FontSize',18,'FontWeight','bold');
            ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
            rmse1=rms(track_selected(:,3)-(interp1(tg_lu(:,1),tg_lu(:,3),track_selected(:,1))),'omitnan');
            rmse2=rms(track_selected(:,3)-(interp1(tg_selected(:,1),tg_selected(:,3),track_selected(:,1))),'omitnan');
            text(tg_lu(1,1),tg_lu(1,3)+50,strcat('RMSE[cm]=',num2str(rmse1)),'Color','r','FontSize',14)
            text(tg_lu(1,1),tg_lu(1,3)+30,strcat('RMSE[cm]=',num2str(rmse2)),'Color',[0.4660 0.6740 0.1880],'FontSize',14)
            yyaxis right
            ylabel('no DAC','Color','k')
            set(gca,'yticklabel',[],'YColor','k')
            hold off
            
            RMSE.nodaclu(H+1)=rmse1;
            RMSE.nodacw(H+1)=rmse2;
            clear rmse1
            clear rmse2
            
            %with dac & oct
            subplot(3,1,3)
            plot(track_selected(:,1),track_selected(:,5),'.','Color',[0.3010 0.7450 0.9330])
            hold on
            plot(track_selected(:,1),smoothdata(track_selected(:,5),'gaussian',50),'.b')
            %with LU
            plot(tg_lu(:,1),tg_lu(:,3),'or')
            plot(tg_lu(:,1),tg_lu(:,3),'--r')
            %withouth LU
            
            plot(tg_selected(:,1),tg_selected(:,3),'o','Color',[0.4660 0.6740 0.1880])
            plot(tg_selected(:,1),tg_selected(:,3),'--','Color',[0.4660 0.6740 0.1880])
            xlabel('Latitude','FontSize',18,'FontWeight','bold');           
            ylabel('DT [cm]','FontSize',18,'FontWeight','bold');
            ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
            rmse1=rms(track_selected(:,5)-(interp1(tg_lu(:,1),tg_lu(:,3),track_selected(:,1))),'omitnan');
            rmse2=rms(track_selected(:,5)-(interp1(tg_selected(:,1),tg_selected(:,3),track_selected(:,1))),'omitnan');
            text(tg_lu(1,1),tg_lu(1,3)+50,strcat('RMSE[cm]=',num2str(rmse1)),'Color','r','FontSize',14)
            text(tg_lu(1,1),tg_lu(1,3)+30,strcat('RMSE[cm]=',num2str(rmse2)),'Color',[0.4660 0.6740 0.1880],'FontSize',14)
            yyaxis right
            ylabel('- OCT','Color','k')
            set(gca,'yticklabel',[],'YColor','k')
            hold off       
            RMSE.octlu(H+1)=rmse1;
            RMSE.octw(H+1)=rmse2;
            clear rmse1
            clear rmse2
            
        end
    end
end
%%
k=id(1);

geoplot(SA_cor.lat(SA_cor.id==k),SA_cor.lon(SA_cor.id==k),'.b')
hold on
geoplot(tg_selected(:,1),tg_selected(:,2),'^k','MarkerFaceColor','k','MarkerSize',10)
text(tg_selected(:,1),tg_selected(:,2)+0.1,num2str(tg_selected(:,4)),'Color','r','FontSize',12)
ax=gca; ax.FontSize=14; ax.FontWeight='Bold'; grid on
set(gca,'fontname','Times New Roman');


%% correlation

CORR=table(); % CORR table

for k=id(1):id(end)
    H = height(CORR);
    
    
    if~isempty(unique(SA_cor.pass(SA_cor.id==k)))
        cycle=unique(SA_cor.cycle(SA_cor.id==k));
        time=SA_cor.date(SA_cor.id==k);
        time=time(1);
        CORR.month(H+1)=month(time);
        CORR.year(H+1)=year((time));
        CORR.cycle(H+1)=cycle;
        
        str=datestr(time,'mmm yyyy');
        
        % Lat SA HDM
        track_selected=[SA_cor.lat(SA_cor.id==k),(SA_cor.ssh(SA_cor.id==k)-SA_cor.N(SA_cor.id==k)+SA_cor.Ell_corr(SA_cor.id==k))*100,...
            (SA_cor.ssh(SA_cor.id==k)-SA_cor.N(SA_cor.id==k)+SA_cor.Ell_corr(SA_cor.id==k)+SA_cor.dac(SA_cor.id==k))*100,...
            (SA_cor.ssh(SA_cor.id==k)-SA_cor.N(SA_cor.id==k)+SA_cor.Ell_corr(SA_cor.id==k)-SA_cor.dac(SA_cor.id==k)+SA_cor.oct(SA_cor.id==k))*100];
        
        if ~isempty(track_selected(~any(ismissing(track_selected(:,3)),2),:))
            
            t=SA_cor.t(SA_cor.id==k);
            t=t(1);
            
            y=2;x=1;
            while tg(n1,y)~=0
                tg_lu(x,:)=[TGinfo.Lat(tg(n1,y)),TGinfo.Lon(tg(n1,y)),interp1(decyear(TGdate),TG_wLU(:,tg(n1,y)),t),tg(n1,y)];
                tg_selected(x,:)=[TGinfo.Lat(tg(n1,y)),TGinfo.Lon(tg(n1,y)),interp1(decyear(TGdate),TG(:,tg(n1,y)),t),tg(n1,y)];
                x=x+1;
                y=y+1;
            end
            
            tg_lu = sortrows(tg_lu,1);
            tg_selected = sortrows(tg_selected,1);
            
            
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
            
            tg_sa=(interp1(tg_selected(:,1),tg_selected(:,3),track_selected(:,1)));
            tg_lusa=(interp1(tg_lu(:,1),tg_lu(:,3),track_selected(:,1)));
            
            [m,~]=find(isnan(tg_sa));
            tg_sa(m,:)=[];
            tg_lusa(m,:)=[];
            track_selected(m,:)=[];
            
            figure(k)
            %without dac
            subplot(3,2,1)
            plot(track_selected(:,2),tg_sa(:,1),'.','Color',[0.4660 0.6740 0.1880])
            cor1=corrcoef(track_selected(:,2),tg_sa(:,1));
            cor1=cor1(1,2);
            title(strcat('No LU cycle#',num2str(cycle),' :',str))
            ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
            ylabel('SA_D_T','FontSize',18,'FontWeight','bold');
            
            subplot(3,2,2)
            plot(track_selected(:,2),tg_lusa(:,1),'.r')
            cor2=corrcoef(track_selected(:,2),tg_lusa(:,1));
            cor2=cor2(1,2);
            ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
            
            CORR.daclu(H+1)=cor2;
            CORR.dacw(H+1)=cor1;
            clear cor1
            clear cor2
            title(strcat('LU cycle#',num2str(cycle),' :',str))
            yyaxis right
            ylabel('+ DAC','Color','k')
            set(gca,'yticklabel',[],'YColor','k')
            
            
            %with dac
            subplot(3,2,3)
            plot(track_selected(:,3),tg_sa(:,1),'.','Color',[0.4660 0.6740 0.1880])
            cor1=corrcoef(track_selected(:,3),tg_sa(:,1));
            cor1=cor1(1,2);
            ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
            ylabel('SA_D_T','FontSize',18,'FontWeight','bold');
            
            subplot(3,2,4)
            plot(track_selected(:,3),tg_lusa(:,1),'.r')
            cor2=corrcoef(track_selected(:,3),tg_lusa(:,1));
            cor2=cor2(1,2);
            ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
            
            CORR.nodaclu(H+1)=cor2;
            CORR.nodacw(H+1)=cor1;
            clear cor1
            clear cor2
            yyaxis right
            ylabel('no DAC','Color','k')
            set(gca,'yticklabel',[],'YColor','k')
            ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
            
            %with dac & OCT
            subplot(3,2,5)
            plot(track_selected(:,4),tg_sa(:,1),'.','Color',[0.4660 0.6740 0.1880])
            cor1=corrcoef(track_selected(:,4),tg_sa(:,1));
            cor1=cor1(1,2);
            ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
            ylabel('SA_D_T','FontSize',18,'FontWeight','bold');
            
            xlabel('TG_D_T','FontSize',18,'FontWeight','bold');
            
            subplot(3,2,6)
            plot(track_selected(:,4),tg_lusa(:,1),'.r')
            cor2=corrcoef(track_selected(:,4),tg_lusa(:,1));
            cor2=cor2(1,2);
            
            CORR.octlu(H+1)=cor2;
            CORR.octacw(H+1)=cor1;
            clear cor1
            clear cor2
            yyaxis right
            ylabel('- OCT','Color','k')
            set(gca,'yticklabel',[],'YColor','k')
            
            xlabel('TG_D_T','FontSize',18,'FontWeight','bold');
            ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
            
            
        end
    end
end


figure(k+1)

subplot(3,1,1)
plot(CORR.cycle,CORR.daclu,'r','LineWidth',3)
hold on
plot(CORR.cycle,CORR.dacw,'Color',[0.4660 0.6740 0.1880],'LineWidth',3)
ylabel('Correlation','FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
ylim([-1 1])
xlim([10 47])

yyaxis right
ylabel('No DAC','Color','k')
set(gca,'yticklabel',[],'YColor','k')
subplot(3,1,2)
plot(CORR.cycle,CORR.nodaclu,'r','LineWidth',3)
hold on
plot(CORR.cycle,CORR.nodacw,'Color',[0.4660 0.6740 0.1880],'LineWidth',3)
ylabel('Correlation','FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
ylim([-1 1])
xlim([10 47])
yyaxis right
ylabel('+ DAC','Color','k')
set(gca,'yticklabel',[],'YColor','k')
subplot(3,1,3)
plot(CORR.cycle,CORR.octlu,'r','LineWidth',3)
hold on
plot(CORR.cycle,CORR.octacw,'Color',[0.4660 0.6740 0.1880],'LineWidth',3)
ylabel('Correlation','FontSize',18,'FontWeight','bold');
xlabel('Cycle','FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
ylim([-1 1])
xlim([10 47])

yyaxis right
ylabel('- OCT','Color','k')
set(gca,'yticklabel',[],'YColor','k')
%% 
DEG=table();
for k=id(1):id(end)
       H = height(DEG);

    
    if~isempty(unique(SA_cor.pass(SA_cor.id==k)))
        cycle=unique(SA_cor.cycle(SA_cor.id==k));
        time=SA_cor.date(SA_cor.id==k);
        time=time(1);
        DEG.month(H+1)=month(time);
        DEG.year(H+1)=year((time));
        DEG.cycle(H+1)=cycle;
        
        str=datestr(time,'mmm yyyy');
        
        % Lat SA HDM
        track_selected=[SA_cor.lat(SA_cor.id==k),...
            (SA_cor.ssh(SA_cor.id==k)-SA_cor.N(SA_cor.id==k)+SA_cor.Ell_corr(SA_cor.id==k)+SA_cor.dac(SA_cor.id==k))*100];
        
        if ~isempty(track_selected(~any(ismissing(track_selected(:,2)),2),:))
            
            t=SA_cor.t(SA_cor.id==k);
            t=t(1);
            
            y=2;x=1;
            while tg(n1,y)~=0
                tg_lu(x,:)=[TGinfo.Lat(tg(n1,y)),TGinfo.Lon(tg(n1,y)),interp1(decyear(TGdate),TG_wLU(:,tg(n1,y)),t),tg(n1,y)];
                tg_selected(x,:)=[TGinfo.Lat(tg(n1,y)),TGinfo.Lon(tg(n1,y)),interp1(decyear(TGdate),TG(:,tg(n1,y)),t),tg(n1,y)];
                x=x+1;
                y=y+1;
            end
            
            tg_lu = sortrows(tg_lu,1);
            tg_selected = sortrows(tg_selected,1);
            
            
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
            
            tg_sa=(interp1(tg_selected(:,1),tg_selected(:,3),track_selected(:,1)));
            tg_lusa=(interp1(tg_lu(:,1),tg_lu(:,3),track_selected(:,1)));
            
            [m,~]=find(isnan(tg_sa));
            tg_sa(m,:)=[];
            tg_lusa(m,:)=[];
            track_selected(m,:)=[];
            
            dh_lu=track_selected(:,2)-tg_lusa;
            
            p1 = polyfit(track_selected(:,1),dh_lu,1);
            f1 = polyval(p1,track_selected(:,1));
            
            
            dh_w=track_selected(:,2)-tg_sa;
            
            
            p2 = polyfit(track_selected(:,1),dh_w,1);
            f2 = polyval(p2,track_selected(:,1));
            
            
            DEG.lu(H+1)=rad2deg(atan(p1(1)));
            DEG.w(H+1)=rad2deg(atan(p2(1)));
            
            figure(k)
            %without dac
%           plot(track_selected(:,1),smoothdata(dh_w,'gaussian',50),'.','Color',[0.4660 0.6740 0.1880])
            plot(track_selected(:,1),dh_w,'.','Color',[.2 1 .2])
            hold on
            plot(track_selected(:,1),f2,'-','Color',[0.4660 0.6740 0.1880],'LineWidth',2.5)
            text(min(track_selected(:,1)+2),max(track_selected(:,2))-20,strcat(sprintf('y = %0.2f x + %0.2f',p2(1),p2(2)),' (tilt=',num2str(rad2deg(atan(p2(1))),2),'\circ)'),'Color',[0.4660 0.6740 0.1880],'FontName','Times New Roman','FontWeight','Bold','FontSize',13)

%           plot(track_selected(:,1),smoothdata(dh_lu,'gaussian',50),'.r')
            plot(track_selected(:,1),dh_lu,'.','Color',[1 .5 .5])
            plot(track_selected(:,1),f1,'-r','LineWidth',2.5)
            text(min(track_selected(:,1)+2),max(track_selected(:,2))-30,strcat(sprintf('y = %0.2f x + %0.2f',p1(1),p1(2)),' (tilt=',num2str(rad2deg(atan(p1(1))),2),'\circ)'),'Color','r','FontName','Times New Roman','FontWeight','Bold','FontSize',13)

            yline(0,'--','Refrence','LineWidth',1.5,'Color',[.3 .3 .3])
            cor1=corrcoef(track_selected(:,2),tg_sa(:,1));
            cor1=cor1(1,2);
            title(strcat('cycle#',num2str(cycle),': ',str))
            ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
            ylabel('\Delta DT [cm]','FontSize',18,'FontWeight','bold');
            pbaspect([2 .5 1])
                        
        end
    end
end


figure(k+2)

subplot(2,2,1:2)
h(1)=plot(DEG.cycle,DEG.lu,'r','LineWidth',3,'displayname',num2str(mean(DEG.lu)));
hold on
h(2)=plot(DEG.cycle,DEG.w,'Color',[0.4660 0.6740 0.1880],'LineWidth',3,'displayname',num2str(mean(DEG.w)));
yline(0,'-','Refrence','LineWidth',1.5)
lgd=legend(h);
title(lgd,'MEAN (\circ)');
ylabel('Tilt \circ','FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
xlabel('Cycle','FontSize',18,'FontWeight','bold');

subplot(2,2,3)
h=histogram(DEG.lu,-87.5:5:87.5);
xlabel('tilt \circ','FontSize',18,'FontWeight','bold');
ylabel('Freq.','FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;

h.FaceColor='r';
subplot(2,2,4)
h=histogram(DEG.w,-87.5:5:87.5);
xlabel('tilt \circ','FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;

h.FaceColor=[0.4660 0.6740 0.1880];