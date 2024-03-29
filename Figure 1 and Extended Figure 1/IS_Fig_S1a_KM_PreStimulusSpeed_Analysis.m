% This code calculate auditory dominance rate according to pre-stimulus
% speed in the keep-moving trial.
% Run this code in the folder where has '01 Stationary Data' and '02 Moving Data' folder

clear variables; close all; clc;

cd('02 Moving Data');

Mat = FindMatFiles();

KMTime = 1.5; % Time epoch to dissociate start-moving or keep-moving trial;
KMThreshold = 3; % Speed threshold for keep-moving trial

for Idx = 1:numel(Mat)
    load(Mat{Idx});

    for j = 1:4
        PreSpeed{Idx,j} = [];
        for ii = 1:size(EventSpeed{4,j},1)
            if nanmean(EventSpeed{4,j}(ii,1:KMTime*EventHz)) > KMThreshold;
                PreSpeed{Idx,j} = cat(1,PreSpeed{Idx,j},nanmean(EventSpeed{4,j}(ii,81:120)));
            end
        end
    end
end

Cut = [3 9 15];

for Idx = 1:size(PreSpeed,1)
    for j = [1 3]
        for k = 1:numel(Cut)
            if k < numel(Cut)
                Ratio{j}(Idx,k) = numel(find(Cut(k) <= PreSpeed{Idx,j} & Cut(k+1) > PreSpeed{Idx,j}))/(numel(find(Cut(k) <= PreSpeed{Idx,j} & Cut(k+1) > PreSpeed{Idx,j})) + numel(find(Cut(k) <= PreSpeed{Idx,j+1} & Cut(k+1) > PreSpeed{Idx,j+1})));
            else
                Ratio{j}(Idx,k) = numel(find(Cut(k) <= PreSpeed{Idx,j}))/(numel(find(Cut(k) <= PreSpeed{Idx,j})) + numel(find(Cut(k) <= PreSpeed{Idx,j+1})));
            end
        end
    end
end

cd ../
%%
xTick1 = [0.9:1:2.9];
xTick2 = [1.1:1:3.1];
Fig = figure('Position',[50 50 110 110]);
hold on;

Data = Ratio{3}.*100;
Avg = nanmean(Data,1); 
for j = 1:size(Data,2)
    Sem(j) = nanstd(Data(:,j),0,1)./sqrt(size(Data(:,j),1)-numnan(Data(:,j)));
end
errorbar(xTick1,Avg,Sem,Sem,'color','k','lineStyle','-','lineWidth',1,'CapSize',3);

Data = Ratio{1}.*100;
Avg = nanmean(Data,1); 
for j = 1:size(Data,2)
    Sem(j) = nanstd(Data(:,j),0,1)./sqrt(size(Data(:,j),1)-numnan(Data(:,j)));
end
errorbar(xTick2,Avg,Sem,Sem,'color','k','lineStyle','-','lineWidth',1,'CapSize',3);

xlim([0.5 3.5])
ylim([0 100])
xticklabels({'3 - 9','9 - 15','15 -'})

set(gca,'TickDir','out','FontName','Arial','FontSize',6);
ylabel('Aud dominance rate (%)','FontName','Arial','FontSize',6);
xlabel('Pre-stimulus speed (cm/s)','FontName','Arial','FontSize',6);

mkdir('Figure'); cd('Figure');
saveas(Fig,'Fig S1a, Dominance and LocomotionSpeed.svg')
cd ../

