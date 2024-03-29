% This code calculates response amplitudes of the incongruent trials
% according to locomotion state (Figure 4j)
% Run this code in the folder where the exmaple data are.

clear variables; close all; clc;

MatFiles = FindMatFiles();
KMTime = 1.5; % Time epoch to sort keep-moving trials
KMThreshold = 3; % speed threshold for sorting keep-moving trials
PreTime = 1;
AmpPre = 0.2;
AmpPost = 1;

nS = 0; nM = 0;
for MatIdx = 1:numel(MatFiles)
    if strcmp(MatFiles{MatIdx}(1:6),'Moving') == 1
        nM = nM + 1;
        MovingMat{nM,1} = MatFiles{MatIdx};
    elseif strcmp(MatFiles{MatIdx}(1:10),'Stationary') == 1
        nS = nS + 1;
        StationaryMat{nS,1} = MatFiles{MatIdx};
    end
end

Target = [2 1; 2 3; 4 1; 4 2; 4 3; 4 4;];

for MatIdx = 1:numel(StationaryMat)
    for Case = 1:size(Target)
        Pre{MatIdx}{Case,1} = [];
        Amp{MatIdx}{Case,1} = [];
        Pre{MatIdx}{Case,3} = [];
        Amp{MatIdx}{Case,3} = [];
        Pre{MatIdx}{Case,2} = [];
        Amp{MatIdx}{Case,2} = [];
    end
    load(StationaryMat{MatIdx});

    for Case = 1:size(Target)
        TempPre1 = []; TempPre2 = [];
        TempAmp1 = []; TempAmp2 = [];
        x = Target(Case,1);
        y = Target(Case,2);
        Avg = nanmean(EventF{x,y},3);
        TempPre1 = nanmean(EventF{x,y}((PSTHPre-PreTime)*ImgHz+1:PSTHPre*ImgHz,:,:),1);
        for k = 1:size(TempPre1,3);
            TempPre2(k,:) = TempPre1(1,:,k);
        end
        Pre{MatIdx}{Case,3} = TempPre2;
        TempAmp1 = nanmean(EventF{x,y}(PSTHPre*ImgHz+1:(PSTHPre+AmpPost)*ImgHz,:,:),1) - nanmean(EventF{x,y}((PSTHPre-AmpPre)*ImgHz:PSTHPre*ImgHz,:,:),1);

        for k = 1:size(TempAmp1,3)
            TempAmp2(k,:) = TempAmp1(1,:,k);
        end
        Temp = TempAmp2(:,:);
        Amp{MatIdx}{Case,3} = Temp;
    end

    load(MovingMat{MatIdx});

    for Case = 1:size(Target)
        x = Target(Case,1);
        y = Target(Case,2);
        KMIdx = []; SMIdx = [];
        for k = 1:size(EventSpeed{x,y},1)
            if nanmean(EventSpeed{x,y}(k,1:KMTime*EventHz),2) >= KMThreshold
                KMIdx = cat(1,KMIdx,k);
            else
                SMIdx = cat(1,SMIdx,k);
            end
        end

        TempPre1 = []; TempPre2 = [];
        TempAmp1 = []; TempAmp2 = [];

        TempPre1 = nanmean(EventF{x,y}((PSTHPre-PreTime)*ImgHz+1:PSTHPre*ImgHz,:,KMIdx),1);
        for k = 1:size(TempPre1,3);
            TempPre2(k,:) = TempPre1(1,:,k);
        end
        Pre{MatIdx}{Case,1} = TempPre2;

        TempAmp1 = nanmean(EventF{x,y}(PSTHPre*ImgHz+1:(PSTHPre+AmpPost)*ImgHz,:,KMIdx),1) - nanmean(EventF{x,y}((PSTHPre-AmpPre)*ImgHz:PSTHPre*ImgHz,:,KMIdx),1);

        for k = 1:size(TempAmp1,3);
            TempAmp2(k,:) = TempAmp1(1,:,k);
        end
        Temp = TempAmp2(:,:);
        Amp{MatIdx}{Case,1} = Temp;


        TempPre1 = []; TempPre2 = [];
        TempAmp1 = []; TempAmp2 = [];

        TempPre1 = nanmean(EventF{x,y}((PSTHPre-PreTime)*ImgHz+1:PSTHPre*ImgHz,:,SMIdx),1);
        for k = 1:size(TempPre1,3);
            TempPre2(k,:) = TempPre1(1,:,k);
        end
        Pre{MatIdx}{Case,2} = TempPre2;

        TempAmp1 = nanmean(EventF{x,y}(PSTHPre*ImgHz+1:(PSTHPre+AmpPost)*ImgHz,:,SMIdx),1) - nanmean(EventF{x,y}((PSTHPre-AmpPre)*ImgHz:PSTHPre*ImgHz,:,SMIdx),1);

        for k = 1:size(TempAmp1,3);
            TempAmp2(k,:) = TempAmp1(1,:,k);
        end
        Temp = TempAmp2(:,:);
        Amp{MatIdx}{Case,2} = Temp;
    end
end

clearvars VAPre AVPre VAAmp AVAmp

for Idx = 1:numel(Pre)
%Idx = 1;

    for i = 3:4
        for j = 1:3
            AVSortedPre{Idx,1}{i-2,j} = [];
            AVSortedPre{Idx,1}{i-2,j} = nanmean(Pre{Idx}{i,j},2);
            if i == 3
                AVSortedPre{Idx,1}{i-2,j} = cat(2,AVSortedPre{Idx}{i-2,j},zeros(numel(AVSortedPre{Idx}{i-2,j}),1)+1);
            elseif i == 4
                AVSortedPre{Idx,1}{i-2,j} = cat(2,AVSortedPre{Idx}{i-2,j},zeros(numel(AVSortedPre{Idx}{i-2,j}),1));
            end
        end
    end

    for j = 1:3
        AVSortedPre{Idx}{1,j} = cat(1,AVSortedPre{Idx}{1,j},AVSortedPre{Idx}{2,j});
    end
    AVSortedPre{Idx}(2,:) = [];


    for i = 5:6
        for j = 1:3
            VASortedPre{Idx,1}{i-4,j} = [];
            VASortedPre{Idx,1}{i-4,j} = nanmean(Pre{Idx}{i,j},2);
            if i == 5
                VASortedPre{Idx,1}{i-4,j} = cat(2,VASortedPre{Idx}{i-4,j},zeros(numel(VASortedPre{Idx}{i-4,j}),1)+1);
            else
                VASortedPre{Idx,1}{i-4,j} = cat(2,VASortedPre{Idx}{i-4,j},zeros(numel(VASortedPre{Idx}{i-4,j}),1));
            end
        end
    end

    for j = 1:3
        VASortedPre{Idx}{1,j} = cat(1,VASortedPre{Idx}{1,j},VASortedPre{Idx}{2,j});
    end
    VASortedPre{Idx}(2,:) = [];
    
    for i = 3:4
        for j = 1:3
            AVSortedAmp{Idx,1}{i-2,j} = [];
            AVSortedAmp{Idx,1}{i-2,j} = nanmean(Amp{Idx}{i,j},2);
            if i == 3
                AVSortedAmp{Idx,1}{i-2,j} = cat(2,AVSortedAmp{Idx}{i-2,j},zeros(numel(AVSortedAmp{Idx}{i-2,j}),1)+1);
            else
                AVSortedAmp{Idx,1}{i-2,j} = cat(2,AVSortedAmp{Idx}{i-2,j},zeros(numel(AVSortedAmp{Idx}{i-2,j}),1));
            end
        end
    end

    for j = 1:3
        AVSortedAmp{Idx}{1,j} = cat(1,AVSortedAmp{Idx}{1,j},AVSortedAmp{Idx}{2,j});
    end
    AVSortedAmp{Idx}(2,:) = [];

    
    for i = 5:6
        for j = 1:3
            VASortedAmp{Idx,1}{i-4,j} = [];
            VASortedAmp{Idx,1}{i-4,j} = nanmean(Amp{Idx}{i,j},2);
            if i == 5
                VASortedAmp{Idx,1}{i-4,j} = cat(2,VASortedAmp{Idx}{i-4,j},zeros(numel(VASortedAmp{Idx}{i-4,j}),1)+1);
            else
                VASortedAmp{Idx,1}{i-4,j} = cat(2,VASortedAmp{Idx}{i-4,j},zeros(numel(VASortedAmp{Idx}{i-4,j}),1));
            end
        end
    end

    for j = 1:3
        VASortedAmp{Idx}{1,j} = cat(1,VASortedAmp{Idx}{1,j},VASortedAmp{Idx}{2,j});
    end
    VASortedAmp{Idx}(2,:) = [];
end

%%
Color = [0 0 0; 0.4 0.4 0.4; 0.85 0.85 0.85];
AvgColor = [1 0 0; 1 0.4 0.4; 1 0.7 0.7];

Data = VASortedAmp;

fig = figure('Position',[50 50 110 110]);
hold on
X = []; Y = [];
for Idx = 1:numel(Data)
    clearvars x y
    for j = 1:numel(Data{Idx})
        x(1,j) = nanmean(Data{Idx}{j}(:,1));
        y(1,j) = nanmean(Data{Idx}{j}(:,2))*100;
    end
    
    plot(x,y,'lineWidth',0.5,'color',[0.9 0.9 0.9]);
    for j = 1:numel(Data{Idx})
        scatter(x(1,j),y(1,j),15,'o','MarkerFaceColor',Color(j,:),'MarkerEdgeColor','none');
        alpha(0.4);
    end
    X = cat(1,X,x);
    Y = cat(1,Y,y);
end

AvgX = nanmean(X,1);
AvgY = nanmean(Y,1);

SemX = nanstd(X,0,1);%./sqrt(size(X,1));
SemY = nanstd(Y,0,1);%./sqrt(size(Y,1));

%plot(AvgX,AvgY,'color','r','lineWidth',1);
for i = 1:3
    errorbar(AvgX(i),AvgY(i),SemY(i),SemY(i),SemX(i),SemX(i),'color',AvgColor(i,:),'lineWidth',1,'CapSize',0);
end

xlim([-0.05 0.5]);
ylim([0 100]);
set(gca,'TickDir','out','FontName','Arial','FontSize',6);
xlabel('Response amplitude, ΔF/F (%)','FontName','Arial','FontSize',6);
ylabel('Aud dominance rate (%)','FontName','Arial','FontSize',6);

mkdir('Figure'); cd('Figure')
saveas(fig,'Fig 4j, VgoAnogo MvmtSorted ResponseAmplitude.svg');
cd ../

x_p(1) = signrank(X(:,1),X(:,2));
x_p(2) = signrank(X(:,2),X(:,3));
x_p(3) = signrank(X(:,1),X(:,3));

y_p(1) = signrank(Y(:,1),Y(:,2));
y_p(2) = signrank(Y(:,2),Y(:,3));
y_p(3) = signrank(Y(:,1),Y(:,3));

Data = AVSortedAmp;

fig = figure('Position',[50 50 110 110]);
hold on
X = []; Y = [];
for Idx = 1:numel(Data)
    clearvars x y
    for j = 1:numel(Data{Idx})
        x(1,j) = nanmean(Data{Idx}{j}(:,1));
        y(1,j) = nanmean(Data{Idx}{j}(:,2))*100;
    end
    
    plot(x,y,'lineWidth',0.5,'color',[0.9 0.9 0.9]);
    for j = 1:numel(Data{Idx})
        scatter(x(1,j),y(1,j),15,'o','MarkerFaceColor',Color(j,:),'MarkerEdgeColor','none');
        alpha(0.4);
    end
    X = cat(1,X,x);
    Y = cat(1,Y,y);
end

AvgX = nanmean(X,1);
AvgY = nanmean(Y,1);

SemX = nanstd(X,0,1);%./sqrt(size(X,1));
SemY = nanstd(Y,0,1);%./sqrt(size(Y,1));

%plot(AvgX,AvgY,'color','r','lineWidth',1);
for i = 1:3
    errorbar(AvgX(i),AvgY(i),SemY(i),SemY(i),SemX(i),SemX(i),'color',AvgColor(i,:),'lineWidth',1,'CapSize',0);
end

xlim([-0.05 0.5]);
ylim([0 100]);
set(gca,'TickDir','out','FontName','Arial','FontSize',6);
xlabel('Response amplitude, ΔF/F (%)','FontName','Arial','FontSize',6);
ylabel('Aud dominance rate (%)','FontName','Arial','FontSize',6);

mkdir('Figure'); cd('Figure')
saveas(fig,'Fig 4j, AgoVnogo MvmtSorted ResponseAmplitude.svg');
cd ../