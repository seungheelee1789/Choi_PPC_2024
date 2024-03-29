%This code makes figure of extended data figure 1d.
%Run this code at folder that has "Pupil Data" folder.
%in the "PupilData.mat' file, there are trial-averaged normalized pupil size data 

clear variables; close all; clc;

cd('03 Pupil Data')

load('PupilData.mat');

cd ../
%Color = {[200/255 200/255 255/255],[255/255 200/255 200/255],[200/255 200/255 255/255],[255/255 200/255 200/255]};

S_Avg = mean(S_Data,1); S_Sem = std(S_Data,0,1)/sqrt(size(S_Data,1));
M_Avg = mean(M_Data,1); M_Sem = std(M_Data,0,1)/sqrt(size(M_Data,1));

XTick = [1 2 3.5 4.5];
XTick2 = XTick + 5;

cd ../

PupilFig = figure('Position',[0 0 150 110]);
hold on

for i = 1:numel(XTick)
    bar(XTick(i),S_Avg(i),'lineWidth',1,'EdgeColor','k','FaceColor','none');
end

for i = 1:size(S_Data,1)
    plot(XTick(1:2),S_Data(i,1:2),'Color',[0.5 0.5 0.5],'lineWidth',0.5);
    plot(XTick(3:4),S_Data(i,3:4),'Color',[0.5 0.5 0.5],'lineWidth',0.5);
end

for i = 1:numel(XTick)
    errorbar(XTick(i),S_Avg(i),S_Sem(i),'CapSize',3,'Color','k','lineWidth',0.75);
end

for i = 1:numel(XTick)
    bar(XTick2(i),M_Avg(i),'lineWidth',1,'EdgeColor','k','FaceColor','none','lineStyle',':');
end

for i = 1:size(M_Data,1)
    plot(XTick2(1:2),M_Data(i,1:2),'Color',[0.5 0.5 0.5],'lineWidth',0.5);
    plot(XTick2(3:4),M_Data(i,3:4),'Color',[0.5 0.5 0.5],'lineWidth',0.5);
end

for i = 1:numel(XTick2)
    errorbar(XTick2(i),M_Avg(i),M_Sem(i),'CapSize',3,'Color','k','lineWidth',0.75);
end

xlim([0 10.5]);
xticks(sort([XTick XTick2],'ascend'));

ylim([0 1]);
yticks(0:0.2:1);
xticklabels({'A','V','A','V','A','V','A','V'});
set(gca,'TickDir','out','FontName','Arial','FontSize',6);
xlabel('ANGVG    AGVNG    ANGVG    AGVNG','FontName','Arial','FontSize',6);
ylabel('Normalized pupil size','FontName','Arial','FontSize',6);

mkdir('Figure'); cd('Figure');
saveas(PupilFig,'Fig S1d, Pre_Stim pupil size.svg');
cd ../
