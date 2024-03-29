% Run this code at the folder where has '01 Stationary Data' and '02 Moving Data' folder

clear variables; close all; clc;

cd('01 Stationary Data');

Mat = FindMatFiles();

for Idx = 1:numel(Mat);

    load(Mat{Idx});
    for i = 1:4
        Stationary(Idx,i) = (TrialNumber(i,1) + TrialNumber(i,3))/sum(TrialNumber(i,:));
    end
    ADomRate{1}(Idx,1) = TrialNumber(4,3)./sum(TrialNumber(4,3:4));
    ADomRate{2}(Idx,1) = TrialNumber(4,1)./sum(TrialNumber(4,1:2));
end

cd ../

%%

clearvars -except Stationary ADomRate; close all; clc;

cd('02 Moving Data');

Mat = FindMatFiles();
KMTime = 1.5; % Time epoch to dissociate start-moving or keep-moving trial;
KMThreshold = 3; % Speed threshold for keep-moving trial

for Idx = 1:numel(Mat);

    load(Mat{Idx});
    for i = 1:4
        Moving(Idx,i) = (TrialNumber(i,1) + TrialNumber(i,3))/sum(TrialNumber(i,:));
    end
    for j = 1:4
        nKM(1,j) = 0; nSM(1,j) = 0;
    end
    for j = 1:4
        for i = 1:size(EventSpeed{4,j},1)
            if nanmean(EventSpeed{4,j}(i,1:EventHz*KMTime)) >= KMThreshold
                nKM(1,j) = nKM(1,j) + 1;
            else
                nSM(1,j) = nSM(1,j) + 1;
            end
        end
    end

    ADomRate{1}(Idx,2) = nSM(1,3)./sum(nSM(1,3:4));
    ADomRate{2}(Idx,2) = nSM(1,1)./sum(nSM(1,1:2));
    ADomRate{1}(Idx,3) = nKM(1,3)./sum(nKM(1,3:4));
    ADomRate{2}(Idx,3) = nKM(1,1)./sum(nKM(1,1:2));
end

cd ../
%%

save('BehaviorCorrectRate.mat','Stationary','Moving','ADomRate');

%%

clear variables; close all; clc;

EdgeColor = {[230/255 0 18/255],[35/255 0 1],[1 0 1],'k'};
FaceColor = {[1 150/255 150/255],[150/255 150/255 1],[1 150/255 1],[150/255 150/255 150/255]};
MkColor = {[230/255 0 18/255],[35/255 0 1],[1 0 1],'k'};

load([pwd '/' 'BehaviorCorrectRate.mat']);

Stationary = Stationary.*100;
Moving = Moving.*100;

BarGraphFig = figure('Position',[0 0 120 110]);
hold on

X1 = [2.1 1 3.2 4.3];
X2 = [2.6 1.5 3.7 4.8];

MouseN = size(Stationary,1);

xticks(sort([X1 X2],'ascend'));

for i = 1:3
    if i < 4
        bar(X1(i),mean(Stationary(:,i),1),0.4,'FaceColor','none','EdgeColor',EdgeColor{i},'lineWidth',1);
    else
        bar([X1(i)],[mean(Stationary(:,i),1) mean(1-Stationary(:,i),1)],0.4,'stacked','FaceColor','none','EdgeColor',EdgeColor{i},'lineWidth',1);
    end
        
    for j = 1:numel(Stationary(:,i))
        XX = [];
        for ii = 1:numel(Moving(:,i))
            XX = [XX; X1(i)+(rand*0.2-0.1)];
        end
    end
    scatter(XX,Stationary(:,i),5,'o','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',MkColor{i},'lineWidth',0.1);
    alpha(0.2)
    errorbar(X1(i),mean(Stationary(:,i),1),std(Stationary(:,i),0,1)./sqrt(MouseN),std(Stationary(:,i),0,1)./sqrt(MouseN),'CapSize',3,'LineWidth',0.75,'Color',EdgeColor{i});
end

for i = 1:3
    if i < 4
        bar(X2(i),mean(Moving(:,i),1),0.4,'FaceColor',FaceColor{i},'EdgeColor',EdgeColor{i},'lineWidth',1);
    else
        bar(X2(i),[mean(Moving(:,i),1) mean(1-Moving(:,i),1)],0.4,'stacked','FaceColor',FaceColor{i},'EdgeColor',EdgeColor{i},'lineWidth',1);
    end
    for j = 1:numel(Moving(:,i))
        XX = [];
        for ii = 1:numel(Moving(:,i))
            XX = [XX; X2(i)+(rand*0.2-0.1)];
        end
    end
    scatter(XX,Moving(:,i),5,'o','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',MkColor{i},'lineWidth',0.1);
    alpha(0.2)
    errorbar(X2(i),mean(Moving(:,i),1),std(Moving(:,i),0,1)./sqrt(MouseN),std(Moving(:,i),0,1)./sqrt(MouseN),'CapSize',3,'LineWidth',0.75,'Color',EdgeColor{i});
end

xlim([0.6 4.1]); 
ylim([50 100]);

xticklabels({'S','M','S','M','S','M','S','M'});
set(gca,'TickDir','out','FontName','Arial','FontSize',6);

xlabel('Types','FontName','Arial','FontSize',6);
ylabel('Correct rate (%)','FontName','Arial','FontSize',6);

mkdir('Figure'); cd('Figure');
saveas(BarGraphFig,'Fig 1f, Treadmill Behavior Result - A,V,Con.svg');
cd ../

InconBarGraphFig = figure('Position',[0 0 60 110]);
hold on

X1 = [2 4 3 1];
X2 = [2.4 4.4 3.4 1.3];

MouseN = size(Stationary,1);

xticks(sort([X1 X2],'ascend'));

for i = 4
    if i < 4
        bar(X1(i),mean(Stationary(:,i),1),0.4,'FaceColor','none','EdgeColor',EdgeColor{i},'lineWidth',1);
    else
        bar([X1(i)],[mean(Stationary(:,i),1) mean(1-Stationary(:,i),1)],0.25,'stacked','FaceColor','none','EdgeColor',EdgeColor{i},'lineWidth',1);
    end
        
    for j = 1:numel(Stationary(:,i))
        XX = [];
        for ii = 1:numel(Moving(:,i))
            XX = [XX; X1(i)+(rand*0.2-0.1)];
        end
    end
    scatter(XX,Stationary(:,i),5,'o','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',MkColor{i},'lineWidth',0.1);
    alpha(0.2)
    errorbar(X1(i),mean(Stationary(:,i),1),std(Stationary(:,i),0,1)./sqrt(MouseN),std(Stationary(:,i),0,1)./sqrt(MouseN),'CapSize',3,'LineWidth',0.75,'Color',EdgeColor{i});
end


for i = 4
    if i < 4
        bar(X2(i),mean(Moving(:,i),1),0.4,'FaceColor',FaceColor{i},'EdgeColor',EdgeColor{i},'lineWidth',1);
    else
        bar(X2(i),[mean(Moving(:,i),1) mean(1-Moving(:,i),1)],0.25,'stacked','FaceColor',FaceColor{i},'EdgeColor',EdgeColor{i},'lineWidth',1);
    end
    for j = 1:numel(Moving(:,i))
        XX = [];
        for ii = 1:numel(Moving(:,i))
            XX = [XX; X2(i)+(rand*0.2-0.1)];
        end
    end
    scatter(XX,Moving(:,i),5,'o','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',MkColor{i},'lineWidth',0.1);
    alpha(0.2)
    errorbar(X2(i),mean(Moving(:,i),1),std(Moving(:,i),0,1)./sqrt(MouseN),std(Moving(:,i),0,1)./sqrt(MouseN),'CapSize',3,'LineWidth',0.75,'Color',EdgeColor{i});
end

xlim([0.6 1.7]); 
ylim([0 100]);

xticklabels({'S','M','S','M','S','M','S','M'});
set(gca,'TickDir','out','FontName','Arial','FontSize',6);

xlabel('Types','FontName','Arial','FontSize',6);
ylabel('Aud dominance rate (%)','FontName','Arial','FontSize',6);

mkdir('Figure'); cd('Figure');
saveas(InconBarGraphFig,'Fig 1g, Treadmill Behavior Result - Incon.svg');
cd ../

%%
%this part is for two-way ANOVA test
% 
% Name = {'Vis','Aud','Con'};
% 
% Data = [];
% DataLabelIdx = 0;
% for j = 1:3
%     for ii = 1:size(Stationary,1)
%         DataLabelIdx = DataLabelIdx + 1;
%         Data = cat(1,Data,Stationary(ii,j).*100);
%         g1{DataLabelIdx,1} = 'S';
%         g2{DataLabelIdx,1} = Name{j};
%     end
% end
% 
% 
% for j = 1:3
%     for ii = 1:size(Moving,1)
%         DataLabelIdx = DataLabelIdx + 1;
%         Data = cat(1,Data,Moving(ii,j).*100);
%         g1{DataLabelIdx,1} = 'M';
%         g2{DataLabelIdx,1} = Name{j};
%     end
% end
% 
% [p,origtbl,stats] = anovan(Data,{g1 g2},"Model","interaction", ...
%     "Varnames",["Loco","Stim"]);
% 
% [results,~,~,gnames] = multcompare(stats,"Dimension",[1 2]);
% 
% tbl = array2table(results,"VariableNames", ...
%     ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
% tbl.("Group A")=gnames(tbl.("Group A"));
% tbl.("Group B")=gnames(tbl.("Group B"))
% 
% tbl = tbl(:,[1 2 6]);

%%

clearvars p;
fig = figure('Position',[50 50 130 110]);
hold on;

Data = ADomRate{1}.*100;

p(1,1) = signrank(Data(:,1),Data(:,2));
p(1,2) = signrank(Data(:,2),Data(:,3));

Avg = nanmean(Data,1);
Sem = nanstd(Data,0,1)./sqrt(size(Data,1));

xTicks1 = [1 2 3];

bar(xTicks1,Avg,'FaceColor','none','EdgeColor','k','LineWidth',1);
for ii = 1:size(Data,1)
    plot(xTicks1,Data(ii,:),'lineWidth',0.5,'color',[180 180 180]./255);
end
errorbar(xTicks1,Avg,Sem,Sem,'lineWidth',0.75,'color','k');

Data = ADomRate{2}.*100;

p(2,1) = signrank(Data(:,1),Data(:,2));
p(2,2) = signrank(Data(:,2),Data(:,3));

Avg = nanmean(Data,1);
Sem = nanstd(Data,0,1)./sqrt(size(Data,1));

xTicks2 = [4.5 5.5 6.5];

bar(xTicks2,Avg,'FaceColor','none','EdgeColor','k','LineWidth',1);
for ii = 1:size(Data,1)
    plot(xTicks2,Data(ii,:),'lineWidth',0.5,'color',[180 180 180]./255);
end
errorbar(xTicks2,Avg,Sem,Sem,'lineWidth',0.75,'color','k');

xlim([0.2 7.3]); 
ylim([0 100]);

xticks([1 2 3 4.5 5.5 6.5]);

xticklabels({'S','SM','KM','S','SM','KM'});
set(gca,'TickDir','out','FontName','Arial','FontSize',6);

xlabel('Types','FontName','Arial','FontSize',6);
ylabel('Aud dominance rate (%)','FontName','Arial','FontSize',6);

mkdir('Figure'); cd('Figure');
saveas(fig,'Fig 1h, Treadmill Auditory Dominance Rate - LocomotionStates.svg');
cd ../