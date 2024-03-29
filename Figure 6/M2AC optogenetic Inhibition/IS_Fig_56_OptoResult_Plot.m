% This code makes figures of optogenetic inhibition during the multisensory behavior task.
% Choose correct file for analysis in line 7-9.

clear variables; close all; clc;

% choose correct file to load
%load(["ACPPC_Inhibition_Opto_TrialNumber.mat"]);
%load(["ACSTR_Inhibition_Opto_TrialNumber.mat"]);
load("M2AC_Inhibition_Opto_TrialNumber.mat");

for i = 1:numel(Main_NoLaserTrial);
    for ii = 1:4
        MainHitRate{1}(i,ii) = Main_NoLaserTrial{i}(ii,1)/sum(Main_NoLaserTrial{i}(ii,1:2))*100;
        MainFARate{1}(i,ii) = Main_NoLaserTrial{i}(ii,4)/sum(Main_NoLaserTrial{i}(ii,3:4))*100;
        MainCR{1}(i,ii) = (Main_NoLaserTrial{i}(ii,1)+Main_NoLaserTrial{i}(ii,3))/sum(Main_NoLaserTrial{i}(ii,1:4))*100;

        MainHitRate{2}(i,ii) = Main_LaserTrial{i}(ii,1)/sum(Main_LaserTrial{i}(ii,1:2))*100;
        MainFARate{2}(i,ii) = Main_LaserTrial{i}(ii,4)/sum(Main_LaserTrial{i}(ii,3:4))*100;
        MainCR{2}(i,ii) = (Main_LaserTrial{i}(ii,1)+Main_LaserTrial{i}(ii,3))/sum(Main_LaserTrial{i}(ii,1:4))*100;

        CtrlHitRate{1}(i,ii) = Ctrl_NoLaserTrial{i}(ii,1)/sum(Ctrl_NoLaserTrial{i}(ii,1:2))*100;
        CtrlFARate{1}(i,ii) = Ctrl_NoLaserTrial{i}(ii,4)/sum(Ctrl_NoLaserTrial{i}(ii,3:4))*100;
        CtrlCR{1}(i,ii) = (Ctrl_NoLaserTrial{i}(ii,1)+Ctrl_NoLaserTrial{i}(ii,3))/sum(Ctrl_NoLaserTrial{i}(ii,1:4))*100;

        CtrlHitRate{2}(i,ii) = Ctrl_LaserTrial{i}(ii,1)/sum(Ctrl_LaserTrial{i}(ii,1:2))*100;
        CtrlFARate{2}(i,ii) = Ctrl_LaserTrial{i}(ii,4)/sum(Ctrl_LaserTrial{i}(ii,3:4))*100;
        CtrlCR{2}(i,ii) = (Ctrl_LaserTrial{i}(ii,1)+Ctrl_LaserTrial{i}(ii,3))/sum(Ctrl_LaserTrial{i}(ii,1:4))*100;
    end
end

MainStim = MainCR{2};
MainNoStim = MainCR{1};
CtrlStim = CtrlCR{2};
CtrlNoStim = CtrlCR{1};

MainDelta = MainStim - MainNoStim;
CtrlDelta = CtrlStim - CtrlNoStim;

EdgeColor = {[230/255 0 18/255],[35/255 0 1],[1 0 1],'k'};
FaceColor = {[1 200/255 200/255],[200/255 200/255 1],[1 200/255 1],[200/255 200/255 200/255]};
MkColor = {[1 0.5 0.5],[0.5 0.5 1],[1 0.5 1],[0.5 0.5 0.5]};
FColor = {[150/255 225/255 150/255],[150/255 225/255 150/255]};

AVCBarGraphFig = figure('Position',[0 0 70 110]);
hold on

X1 = [2 1 3];
X2 = [2.4 1.4 3.4];

MouseN = size(MainDelta,1);
CtrlMouseN = size(CtrlDelta,1);

xticks(sort(mean([X1;X2],1),'ascend'));


for i = 1:numel(X1)
    if i < 4
        bar(X1(i),mean(MainDelta(:,i),1),0.4,'FaceColor',FColor{1},'EdgeColor',EdgeColor{i},'LineWidth',1);
    else
        bar([X1(i)],[mean(MainDelta(:,i),1) mean(1-MainDelta(:,i),1)],0.4,'stacked','FaceColor',FColor{1},'EdgeColor',EdgeColor{i});
    end
    
    if i < 4
        bar(X2(i),mean(CtrlDelta(:,i),1),0.4,'FaceColor','none','EdgeColor',EdgeColor{i},'LineWidth',1);
    else
        bar([X2(i)],[mean(CtrlDelta(:,i),1) mean(1-CtrlDelta(:,i),1)],0.4,'stacked','FaceColor','none','EdgeColor',EdgeColor{i});
    end

    errorbar(X1(i),mean(MainDelta(:,i),1),std(MainDelta(:,i),0,1)./sqrt(MouseN),std(MainDelta(:,i),0,1)./sqrt(MouseN),'CapSize',3,'LineWidth',0.75,'Color',EdgeColor{i});
    errorbar(X2(i),mean(CtrlDelta(:,i),1),std(CtrlDelta(:,i),0,1)./sqrt(CtrlMouseN),std(CtrlDelta(:,i),0,1)./sqrt(CtrlMouseN),'CapSize',3,'LineWidth',0.75,'Color',EdgeColor{i});
end

for i = 1:3
    for j = 1:size(MainDelta,1)
        plot([X1(i)+0.05 X2(i)-0.05],[MainDelta(j,i) CtrlDelta(j,i)],'lineWidth',0.5,'Color',[0.5 0.5 0.5]);
    end
end

xlim([0.6 3.8]); 
ylim([-40 40]);


xticklabels({'Aud' 'Vis' 'Con'});
set(gca,'TickDir','out','FontName','Arial','FontSize',6);

xlabel('Types','FontName','Arial','FontSize',6);
ylabel('ΔCorrect rate (laser on - off, %)','FontName','Arial','FontSize',6);

mkdir('Figure'); cd('Figure');
saveas(AVCBarGraphFig,'01. A,V,Con Correct rate.svg');
cd ../

% signrank(MainDelta(:,1))
% signrank(MainDelta(:,2))
% signrank(MainDelta(:,3))
% 
% signrank(CtrlDelta(:,1))
% signrank(CtrlDelta(:,2))
% signrank(CtrlDelta(:,3))
%

InconBarGraphFig = figure('Position',[0 0 50 110]);
hold on

X2 = [1];
X3 = [1.4];

MouseN = size(MainDelta,1);

xticks(sort(mean([X2 X3]),'ascend'));


for i = 4
    bar([X2(1)],mean(MainDelta(:,i),1),0.3,'stacked','FaceColor',FColor{1},'EdgeColor',EdgeColor{i},'lineWidth',1);
    bar([X3(1)],mean(CtrlDelta(:,i),1),0.3,'stacked','FaceColor','none','EdgeColor',EdgeColor{i},'lineWidth',1);
    errorbar(X2(1),mean(MainDelta(:,i),1),std(MainDelta(:,i),0,1)./sqrt(MouseN),std(MainDelta(:,i),0,1)./sqrt(MouseN),'CapSize',3,'LineWidth',0.75,'Color',EdgeColor{i});
    errorbar(X3(1),mean(CtrlDelta(:,i),1),std(CtrlDelta(:,i),0,1)./sqrt(CtrlMouseN),std(CtrlDelta(:,i),0,1)./sqrt(CtrlMouseN),'CapSize',3,'LineWidth',0.75,'Color',EdgeColor{i});
end

for j = 1:size(MainDelta,1)
    plot([X2(1)+0.05 X3(1)-0.05],[MainDelta(j,4) CtrlDelta(j,4)],'lineWidth',0.5,'Color',[0.5 0.5 0.5]);
end

xlim([0.6 1.8]); 
ylim([-60 60]);

xticklabels({'Incon'});
set(gca,'TickDir','out','FontName','Arial','FontSize',6);

xlabel('Types','FontName','Arial','FontSize',6);
ylabel('ΔAud dominance rate (laser on - off, %)','FontName','Arial','FontSize',6);

mkdir('Figure'); cd('Figure');
saveas(InconBarGraphFig,'02. Incon Correct rate.svg');
cd ../

