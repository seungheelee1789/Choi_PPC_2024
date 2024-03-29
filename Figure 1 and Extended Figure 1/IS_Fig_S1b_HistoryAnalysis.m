%This codes analyze how previous trial affects auditory dominance rate
% Run this code in the folder where has '01 Stationary Data' and '02 Moving Data' folder

clear variables; close all;  clc;

%%
%Collect Stationary Session Data

cd('01 Stationary Data');

Mat = FindMatFiles();

for Idx = 1:numel(Mat)

    load(Mat{Idx})

    UIdx = find(EventSeq(:,1) == 212);
    if UIdx(end) == size(EventSeq,1)
        UIdx = UIdx(1:end-1);
    end
    Target = EventSeq(UIdx+1,:);
    Incon = Target(find(Target(:,1) == 412),:);

    if isempty(Incon) == 0
        ADomRate{1}(Idx,1) = sum(Incon(:,4))./size(Incon,1)*100;
    else
        ADomRate{1}(Idx,1) = nan;
    end

    UIdx = find(EventSeq(:,1) == 111);
    if UIdx(end) == size(EventSeq,1)
        UIdx = UIdx(1:end-1);
    end
    Target = EventSeq(UIdx+1,:);
    Incon = Target(find(Target(:,1) == 412),:);
    
    if isempty(Incon) == 0
        ADomRate{1}(Idx,2) = sum(Incon(:,4))./size(Incon,1)*100;
    else
        ADomRate{1}(Idx,2) = nan;
    end


    UIdx = find(EventSeq(:,1) == 211);
    if UIdx(end) == size(EventSeq,1)
        UIdx = UIdx(1:end-1);
    end
    Target = EventSeq(UIdx+1,:);
    Incon = Target(find(Target(:,1) == 411),:);

    if isempty(Incon) == 0
        ADomRate{2}(Idx,1) = sum(Incon(:,4))./size(Incon,1)*100;
    else
        ADomRate{2}(Idx,1) = nan;
    end

    UIdx = find(EventSeq(:,1) == 112);
    if UIdx(end) == size(EventSeq,1)
        UIdx = UIdx(1:end-1);
    end
    Target = EventSeq(UIdx+1,:);
    Incon = Target(find(Target(:,1) == 411),:);
    
    if isempty(Incon) == 0
        ADomRate{2}(Idx,2) = sum(Incon(:,4))./size(Incon,1)*100;
    else
        ADomRate{2}(Idx,2) = nan;
    end
end

cd ../

%%
%Collect moving session data
clearvars -except ADomRate; close all;  clc;

cd('02 Moving Data');

Mat = FindMatFiles();

for Idx = 1:numel(Mat)

    load(Mat{Idx})

    UIdx = find(EventSeq(:,1) == 212);
    if UIdx(end) == size(EventSeq,1)
        UIdx = UIdx(1:end-1);
    end
    Target = EventSeq(UIdx+1,:);
    Incon = Target(find(Target(:,1) == 412),:);

    if isempty(Incon) == 0
        ADomRate{3}(Idx,1) = sum(Incon(:,4))./size(Incon,1)*100;
    else
        ADomRate{3}(Idx,1) = nan;
    end

    UIdx = find(EventSeq(:,1) == 111);
    if UIdx(end) == size(EventSeq,1)
        UIdx = UIdx(1:end-1);
    end
    Target = EventSeq(UIdx+1,:);
    Incon = Target(find(Target(:,1) == 412),:);
    
    if isempty(Incon) == 0
        ADomRate{3}(Idx,2) = sum(Incon(:,4))./size(Incon,1)*100;
    else
        ADomRate{3}(Idx,2) = nan;
    end


    UIdx = find(EventSeq(:,1) == 211);
    if UIdx(end) == size(EventSeq,1)
        UIdx = UIdx(1:end-1);
    end
    Target = EventSeq(UIdx+1,:);
    Incon = Target(find(Target(:,1) == 411),:);

    if isempty(Incon) == 0
        ADomRate{4}(Idx,1) = sum(Incon(:,4))./size(Incon,1)*100;
    else
        ADomRate{4}(Idx,1) = nan;
    end

    UIdx = find(EventSeq(:,1) == 112);
    if UIdx(end) == size(EventSeq,1)
        UIdx = UIdx(1:end-1);
    end
    Target = EventSeq(UIdx+1,:);
    Incon = Target(find(Target(:,1) == 411),:);
    
    if isempty(Incon) == 0
        ADomRate{4}(Idx,2) = sum(Incon(:,4))./size(Incon,1)*100;
    else
        ADomRate{4}(Idx,2) = nan;
    end
end

cd ../

%%

xTick = [1 2; 3.5 4.5; 6 7; 8.5 9.5];

fig = figure('Position',[50 50 150 110]);
hold on
for j = 1:4
    Data = ADomRate{j};

    p(j) = signrank(Data(:,1),Data(:,2));

    Avg = nanmean(Data,1);
    for jj = 1:2
        Sem(jj) = nanstd(Data(:,jj),0,1)./sqrt(size(Data(:,jj),1)-numnan(Data(:,jj)));
    end

    if j < 3
        bar(xTick(j,1),Avg(1),'lineWidth',1,'EdgeColor','k','FaceColor','none');
        bar(xTick(j,2),Avg(2),'lineWidth',1,'EdgeColor','k','FaceColor','none');
    else
        bar(xTick(j,1),Avg(1),'lineWidth',1,'EdgeColor','k','FaceColor','none','lineStyle',':');
        bar(xTick(j,2),Avg(2),'lineWidth',1,'EdgeColor','k','FaceColor','none','lineStyle',':');
    end

    errorbar(xTick(j,1),Avg(1),Sem(1),Sem(1),'lineWidth',0.75,'CapSize',3,'color','k');
    errorbar(xTick(j,2),Avg(2),Sem(2),Sem(2),'lineWidth',0.75,'CapSize',3,'color','k');

    % for ii = 1:size(Data,1)
    %     plot([xTick(j,1)+0.1 xTick(j,2)-0.1],Data(ii,:),'lineWidth',0.5,'color',[0.5 0.5 0.5]);
    % end

end

xticks([1 2 3.5 4.5 6 7 8.5 9.5]);
xticklabels({'Ang','Vg','Ag','Vng','Ang','Vg','Ag','Vng'});

xlim([0 10.5])
ylim([0 100])
%xticklabels({'3 - 9','9 - 15','15 -'})

set(gca,'TickDir','out','FontName','Arial','FontSize',6);
ylabel('Aud dominance rate (%)','FontName','Arial','FontSize',6);
xlabel('ANGVG    AGVNG   ANGVG    AGVNG','FontName','Arial','FontSize',6);

mkdir('Figure'); cd('Figure');
saveas(fig,'Fig S1b, Trial History.svg')
cd ../

