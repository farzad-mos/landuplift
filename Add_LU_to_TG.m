% add land uplift correction to TG data by NKG2016LU model
% code by VJ --- last Update: 5 Oct 2021
% --------------------------------------------------------

clear
load NKG2016LU_lev.dat
load TGdata.mat %TGinfo TG TGdate

TG552=readtable('552.txt');
TG552.Properties.VariableNames{1} = 'Time';
TG552.Properties.VariableNames{2} = 'obs';
TGID=552;
date=TG552.Time;
TG=TG552.obs;


epoch = 2000;

 
if TGID==552
%Ratan TG only
Lon=20.8950;
Lat=63.9861;
else
Lon=TGinfo.Lon(TGID);
Lat=TGinfo.Lat(TGID);
end



LU = griddata(NKG2016LU_lev(:,2), NKG2016LU_lev(:,1), NKG2016LU_lev(:,3),...
                    Lon,  Lat, "linear")/10; %convert to cm/year
% LU (TGinfo.Country == "Estonia") = 0;

% TGinfo.NKGU2016LU = LU;


DecTGdate = decyear(date);

LUcorrection = (DecTGdate-epoch) * LU';

TG_wLU = TG + LUcorrection;
% save("TGdata.mat","TG_wLU","TG","TGdate","TGinfo")


TG552.LU=TG_wLU;
save("TG552.mat")