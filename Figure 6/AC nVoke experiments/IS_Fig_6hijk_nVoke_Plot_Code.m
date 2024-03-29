% This code plots results of optogenetic manupulation experiments in Figure
% 6h-6k.
% Run this code where the "AC_nVokeData.mat" file is

% This section plots result of ACPPC data
clear variables; close all; clc;

load('AC_nVokeData.mat');

AmpPre = 0.2;
AmpPost = 1;
nPre = PSTHPre*ImgHz;
nPost = PSTHPost*ImgHz;

xAxis = -PSTHPre+1/ImgHz:1/ImgHz:PSTHPost;
yLim = [-0.2 1];

Color = [0 152 243; 0 0 255]./255;
Style{1} = '-'; Style{2} = ':';

Name = {'5kHz','10kHz'};

Merged_FData = ACPPC_Merged_FData;

for i = 1:2
    for j = 1:2
        for k = 1:size(Merged_FData{i,j},3)
            Amp{i,j}(k,:) = nanmean(Merged_FData{i,j}(nPre+1:nPre+AmpPost*ImgHz,:,k),1) - nanmean(Merged_FData{i,j}(nPre-AmpPre*ImgHz:nPre,:,k),1);
        end
        AvgF{i,j} = nanmean(nanmean(Merged_FData{i,j},3),2);
        SemF{i,j} = nanstd(nanmean(Merged_FData{i,j},3),0,2)./sqrt(size(nanmean(Merged_FData{i,j},3),2));

        AvgAmp{i,j} = nanmean(Amp{i,j},1);
    end
end

for i = 1:2
    fig = figure('Position',[50 50 60 110]);
    hold on

    for j = 1:2
        errorshade(xAxis,AvgF{i,j}',AvgF{i,j}'+SemF{i,j}',AvgF{i,j}'-SemF{i,j}',Color(i,:));
        alpha(0.2);
        plot(xAxis,AvgF{i,j},'lineWidth',1,'color',Color(i,:),'lineStyle',Style{j});
    end

    xlim([-1 2]);
    ylim([-0.2 1]);
    yticks([yLim(1):0.2:yLim(2)]);

    line([0 0],yLim,'color',[160 160 160]./255,'lineWidth',0.5);
    line([1 1],yLim,'color',[160 160 160]./255,'lineWidth',0.5);

    set(gca,'TickDir','out','FontName','Arial','FontSize',6);
    xlabel('Time (s)','FontName','Arial','FontSize',6);
    ylabel('ΔF/F (%)','FontName','Arial','FontSize',6);

    mkdir('Figure'); cd('Figure');
    saveas(fig,['Fig 6h, ACPPC ' Name{i} ' PSTH.svg']);
    cd ../
end
%
Order = [3.25 4.25 1 2];

for i = 1:2
    for j = 1:2
        Avg(i,j) = nanmean(AvgAmp{i,j});
        Sem(i,j) = nanstd(AvgAmp{i,j})./sqrt(numel(AvgAmp{i,j}));
    end
end

fig = figure('Position',[50 50 80 110]);
hold on;
for i = 1:2
    for j = 1:2
        bar(Order((i-1)*2+j),Avg(i,j),0.65,'FaceColor','none','EdgeColor',Color(i,:),'lineWidth',1,'lineStyle',Style{j});
        errorbar(Order((i-1)*2+j),Avg(i,j),Sem(i,j),Sem(i,j),'lineWidth',0.75,'color',Color(i,:),'CapSize',3);
    end
end

xlim([0.2 5.05]);
ylim([0 0.8]);

xticks(sort(Order,'ascend'));
xticklabels({'Off','On','Off','On'});

set(gca,'TickDir','out','FontName','Arial','FontSize',6);
xlabel('Type','FontName','Arial','FontSize',6);
ylabel('Response amplitude, ΔF/F (%)','FontName','Arial','FontSize',6);

mkdir('Figure'); cd('Figure');
saveas(fig,'Fig 6i, ACPPC Response amplitude.svg');
cd ../

ACPPC_5kHz_PValue = signrank(AvgAmp{1,1},AvgAmp{1,2})
ACPPC_10kHz_PValue = signrank(AvgAmp{2,1},AvgAmp{2,2})
%%
% This section plots result of ACSTR data
clear variables; close all; 

load('AC_nVokeData.mat');

AmpPre = 0.2;
AmpPost = 1;
nPre = PSTHPre*ImgHz;
nPost = PSTHPost*ImgHz;

xAxis = -PSTHPre+1/ImgHz:1/ImgHz:PSTHPost;
yLim = [-0.2 1];

Color = [0 152 243; 0 0 255]./255;
Style{1} = '-'; Style{2} = ':';

Name = {'5kHz','10kHz'};

Merged_FData = ACSTR_Merged_FData;

for i = 1:2
    for j = 1:2
        for k = 1:size(Merged_FData{i,j},3)
            Amp{i,j}(k,:) = nanmean(Merged_FData{i,j}(nPre+1:nPre+AmpPost*ImgHz,:,k),1) - nanmean(Merged_FData{i,j}(nPre-AmpPre*ImgHz:nPre,:,k),1);
        end
        AvgF{i,j} = nanmean(nanmean(Merged_FData{i,j},3),2);
        SemF{i,j} = nanstd(nanmean(Merged_FData{i,j},3),0,2)./sqrt(size(nanmean(Merged_FData{i,j},3),2));

        AvgAmp{i,j} = nanmean(Amp{i,j},1);
    end
end

for i = 1:2
    fig = figure('Position',[50 50 60 110]);
    hold on

    for j = 1:2
        errorshade(xAxis,AvgF{i,j}',AvgF{i,j}'+SemF{i,j}',AvgF{i,j}'-SemF{i,j}',Color(i,:));
        alpha(0.2);
        plot(xAxis,AvgF{i,j},'lineWidth',1,'color',Color(i,:),'lineStyle',Style{j});
    end

    xlim([-1 2]);
    ylim([-0.2 1]);
    yticks([yLim(1):0.2:yLim(2)]);

    line([0 0],yLim,'color',[160 160 160]./255,'lineWidth',0.5);
    line([1 1],yLim,'color',[160 160 160]./255,'lineWidth',0.5);

    set(gca,'TickDir','out','FontName','Arial','FontSize',6);
    xlabel('Time (s)','FontName','Arial','FontSize',6);
    ylabel('ΔF/F (%)','FontName','Arial','FontSize',6);

    mkdir('Figure'); cd('Figure');
    saveas(fig,['Fig 6j, ACSTR ' Name{i} ' PSTH.svg']);
    cd ../
end
%
Order = [3.25 4.25 1 2];

for i = 1:2
    for j = 1:2
        Avg(i,j) = nanmean(AvgAmp{i,j});
        Sem(i,j) = nanstd(AvgAmp{i,j})./sqrt(numel(AvgAmp{i,j}));
    end
end

fig = figure('Position',[50 50 80 110]);
hold on;
for i = 1:2
    for j = 1:2
        bar(Order((i-1)*2+j),Avg(i,j),0.65,'FaceColor','none','EdgeColor',Color(i,:),'lineWidth',1,'lineStyle',Style{j});
        errorbar(Order((i-1)*2+j),Avg(i,j),Sem(i,j),Sem(i,j),'lineWidth',0.75,'color',Color(i,:),'CapSize',3);
    end
end

xlim([0.2 5.05]);
ylim([0 0.8]);

xticks(sort(Order,'ascend'));
xticklabels({'Off','On','Off','On'});

set(gca,'TickDir','out','FontName','Arial','FontSize',6);
xlabel('Type','FontName','Arial','FontSize',6);
ylabel('Response amplitude, ΔF/F (%)','FontName','Arial','FontSize',6);

mkdir('Figure'); cd('Figure');
saveas(fig,'Fig 6k, ACSTR Response amplitude.svg');
cd ../

ACSTR_5kHz_PValue = signrank(AvgAmp{1,1},AvgAmp{1,2})
ACPPC_10kHz_PValue = signrank(AvgAmp{2,1},AvgAmp{2,2})

