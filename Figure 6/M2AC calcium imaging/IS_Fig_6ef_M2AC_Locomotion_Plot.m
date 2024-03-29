% This code calculate Pearson correlation value between M2AC activity and
% locomotion speed, and make plots of figure 6e and 6f.
% Run this code at where the example data are.
%%
clear variables; close all; clc;

LongLocomotionTime = 5; %if locomotion bout is longer than this value (second), calculate response

Mat = FindMatFiles();

EarlyAvg = []; LateAvg = []; EarlySem = []; LateSem = [];
 OnsetR = [];  OffsetR = []; MergedLongF = [];
for Session = 1:numel(Mat)
    load(Mat{Session});
    
    MergedF = []; MergedSpeed = [];
    for i = 1:numel(LocomotionOnsetSpeed)
        Data = LocomotionOnsetSpeed{i};
        MergedSpeed = cat(2,MergedSpeed,nanmean(reshape([Data(:); nan(mod(-numel(Data),2),1)],2,[]))); % downsampling
        MergedF = cat(1,MergedF,LocomotionOnsetF{i});
    end
    for jj = 1:size(MergedF,2)
        OnsetR = cat(1,OnsetR,corr(MergedSpeed',MergedF(:,jj)));
    end

    MergedF = []; MergedSpeed = [];
    for i = 1:numel(LocomotionOffsetSpeed)
        Data = LocomotionOffsetSpeed{i};
        MergedSpeed = cat(2,MergedSpeed,nanmean(reshape([Data(:); nan(mod(-numel(Data),2),1)],2,[]))); % downsampling
        MergedF = cat(1,MergedF,LocomotionOffsetF{i});
    end
    for jj = 1:size(MergedF,2)
        OffsetR = cat(1,OffsetR,corr(MergedSpeed',MergedF(:,jj)));
    end
    
    clearvars Idx BaseLine Amp1 Amp2
    Idx = find(Duration > LongLocomotionTime);
    LongF = [];
    for i = 1:numel(Idx)
        BaseLine = nanmean(LocomotionBoutF{Idx(i)}(1:1*ImgHz,:),1);
        Amp1(i,:) = nanmean(LocomotionBoutF{Idx(i)}((PreTime)*ImgHz+1:(PreTime+1.5)*ImgHz,:),1) - BaseLine;
        Amp2(i,:) = nanmean(LocomotionBoutF{Idx(i)}((PreTime+3)*ImgHz+1:end-ImgHz*3,:),1) - BaseLine;
        LongF = cat(3,LongF,LocomotionBoutF{Idx(i)}(1:7*ImgHz+1,:));
    end
    EarlyAvg = cat(2,EarlyAvg,nanmean(Amp1,1));
    LateAvg = cat(2,LateAvg,nanmean(Amp2,1));
    MergedLongF = cat(2,MergedLongF,nanmean(LongF,3));
end

signrank(EarlyAvg,LateAvg)

mkdir('Figure'); cd('Figure');

Bin = -1:0.1:1;

fig = figure('Position', [50 50 80 110]); hold on;

% histogram(IncR,Bin,'FaceColor',IncColor,'EdgeColor','k');
% alpha(0.5);
histogram(OnsetR,Bin,'FaceColor','none','EdgeColor','k');
%alpha(0.5);
% scatter(nanmean(IncR),10,'v','MarkerFaceColor',IncColor,'MarkerEdgeColor','none');
scatter(nanmean(OnsetR),15,'v','MarkerFaceColor','r','MarkerEdgeColor','k');

xlim([-1 1]); ylim([0 20]);
%ylim(YAxis);
set(gca,'TickDir','out','FontName','Arial','FontSize',6);
xlabel('Corr. coef, r','FontName','Arial','FontSize',6);
ylabel('Cell number','FontName','Arial','FontSize',6);

saveas(fig,['Fig 6d, Locomotion Onset Pearson Correlation Histogram.svg']);

%%

xAxis = -2:1/ImgHz:5;
TraceFig = figure('Position', [50 50 80 80]); hold on
AvgF = nanmean(MergedLongF,2); SemF = nanstd(MergedLongF,0,2)./sqrt(size(MergedLongF,2));
errorshade(xAxis,AvgF',AvgF'+SemF',AvgF'-SemF',[0 0.75 0]);
alpha(0.2);
plot(xAxis,AvgF,'color',[0 0.75 0]);
xticks([-2:5]);
set(gca,'TickDir','out','FontName','Arial','FontSize',6);

ylim([-1 4]); ylabel('ΔF/F(%)','FontName','Arial','FontSize',6); xlabel('Time (s)','FontName','Arial','FontSize',6);
line([0 0],[-1 4],'lineWidth',0.75,'color',[0.5 0.5 0.5]);
line([1.5 1.5],[-1 4],'lineWidth',0.75,'color',[0.5 0.5 0.5]);
line([3 3],[-1 4],'lineWidth',0.75,'color',[0.5 0.5 0.5]);

saveas(TraceFig,'Fig 6f, Locomotion Onset PSTH.svg');

%%

fig = figure('Position', [50 50 80 80]); hold on;

Data = EarlyAvg;
AvgPlot = nanmean(Data); SemPlot = nanstd(Data)/sqrt(numel(Data));
bar(1,AvgPlot,0.5,'lineWidth',1,'EdgeColor','k','FaceColor','none');
errorbar(1,AvgPlot,SemPlot,SemPlot,'CapSize',3,'linewidth',0.75,'color','k');

Data = LateAvg;
AvgPlot = nanmean(Data); SemPlot = nanstd(Data)/sqrt(numel(Data));
bar(2,AvgPlot,0.5,'lineWidth',1,'EdgeColor','k','FaceColor','none');
errorbar(2,AvgPlot,SemPlot,SemPlot,'CapSize',3,'linewidth',0.75,'color','k');

xlim([0.4 2.6]); ylim([0 5]);
%ylim(YAxis);
set(gca,'TickDir','out','FontName','Arial','FontSize',6);
xticks([1 2]);
xticklabels({'Early','Late'});
xlabel('Time from locomotion onset (s)','FontName','Arial','FontSize',6);
ylabel('ΔM2AC activity, ΔF/F(%)','FontName','Arial','FontSize',6);

saveas(fig,['Fig 6f, Locomotion Onset Activity Change.svg']);

cd ../
