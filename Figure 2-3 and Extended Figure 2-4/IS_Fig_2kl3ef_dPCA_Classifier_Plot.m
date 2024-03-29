% This code make figures about classifier in Figure 2k, 2l, 3e, 3f and
% supplementary figure 3c, 3d
% You should have run "IS_Run_dPCA.m" before this
% Run this code where the "dPCA Result.mat" file is

clear variables; close all; clc;

for Type = 1:1 %1 for Stimulus axes, 2 for Decision axes, 3 for Independent axes, 4 for Interaction axes

    Session = 0;

    load('dPCA Result.mat');

    nShuffle = 500; % number of iteration for shuffling of trial label

    for S_Idx = 1:size(AvgScore,1)
        if isempty(AvgScore{S_Idx}) == 0
            Session = Session + 1;

            disp(['Analyzing ' num2str(Session) 'th data']);
            Data = AvgScore{S_Idx,Type};

            PC = 1:3 %remove this line to use all PCs

            if exist('PC')
                if size(Data{1},2) < max(PC)
                    PC = 1:size(Data{1},2);
                end
            else
                PC = 1:size(Data{1},2);
            end

            Target = [1 4 5 8]; % 1 = Vis hit; 4 = Vis CR; 5 = Aud hit; 8 = Aud CR;
            CatUni = [];
            CatPre = [];
            for Idx = 1:4
                AvgSingle{Idx} = Data{Target(Idx)}(:,PC);
                CatUni = cat(1,CatUni,AvgSingle{Idx});
                CatPre = cat(1,CatPre,AvgSingle{Idx}(1:ImgHz*Pre,:));
            end

            RMS = 1;
            PreRMS = 1;

            AvgCon{1} = Data{9}(:,PC); AvgCon{2} = Data{12}(:,PC);

            for Idx = 1:4
                AvgIncon{Idx} = Data{Idx+12}(:,PC);
            end

            Data = UniTrialScore{S_Idx,Type};
            CatDecodeResult = [];
            for Idx = 1:3
                CatGNGDecodeResult{Idx} = [];
            end

            for Idx = 1:4
                TrialData{Idx} = Data{Target(Idx)}(:,PC,:);
                TrialNum(Idx) = size(TrialData{Idx},3);
                for k = 1:size(TrialData{Idx},3)
                    for i = 1:size(TrialData{Idx},1)
                        D = [];
                        for j = 1:numel(AvgSingle)
                            D(j) = norm(TrialData{Idx}(i,:,k)-AvgSingle{j}(i,:));
                            Uni_DtoUniAvg{S_Idx,Idx}(i,j,k) = D(j)./RMS;
                        end
                        D(5) = norm(TrialData{Idx}(i,:,k)-nanmean(CatPre,1));
                        [~, SingleClass{Session}{Idx}(i,k)] = min(D);

                        D = [];
                        if Idx < 3
                            for j = 1:2
                                D(j) = norm(TrialData{Idx}(i,:,k)-AvgSingle{j}(i,:));
                            end
                            [~,GNG_Class{Session}{Idx}(i,k)] = min(D);
                        elseif Idx > 2
                            D(1:2) = nan;
                            for j = 3:4
                                D(j) = norm(TrialData{Idx}(i,:,k)-AvgSingle{j}(i,:));
                            end
                            [~,GNG_Class{Session}{Idx}(i,k)] = min(D);
                        end
                    end
                    CatDecodeResult = cat(2,CatDecodeResult,SingleClass{Session}{Idx}(:,k));
                    if Idx < 3
                        CatGNGDecodeResult{1} = cat(2,CatGNGDecodeResult{1},GNG_Class{Session}{Idx}(:,k));
                    elseif Idx > 2
                        CatGNGDecodeResult{2} = cat(2,CatGNGDecodeResult{2},GNG_Class{Session}{Idx}(:,k));
                    end
                end
            end

            for N = 1:nShuffle
                if rem(N,100) == 0
                    disp(['Currently at ' num2str(N) 'th shuffling of single trial']);
                end
                Shuffle = randsample(sum(TrialNum),sum(TrialNum),'false');
                for Idx = 1:4
                    if Idx == 1
                        Trials = Shuffle(1:TrialNum(Idx));
                    else
                        Trials = Shuffle(sum(TrialNum(1:Idx-1))+1:sum(TrialNum(1:Idx)));
                    end
                    for ii = 1:size(CatDecodeResult,1)
                        ShuffleCR{Session}(ii,Idx,N) = numel(find(CatDecodeResult(ii,Trials) == Idx)) / TrialNum(Idx);
                    end
                end
            end

            Data = MultiTrialScore{S_Idx,Type};
            Target = [1 4];
            for Idx = 1:numel(Target)
                TrialData{Idx} = Data{Target(Idx)}(:,PC,:);
                TrialNum(Idx+4) = size(TrialData{Idx},3);
                for k = 1:size(TrialData{Idx},3)
                    for i = 1:size(TrialData{Idx},1)
                        D = []; D(1:4) = nan;
                        for j = 1:numel(AvgCon)
                            D(j+4) = norm(TrialData{Idx}(i,:,k)-AvgCon{j}(i,:));
                        end
                        [~, GNG_Class{Session}{Idx+4}(i,k)] = min(D);
                    end
                    CatGNGDecodeResult{3} = cat(2,CatGNGDecodeResult{3},GNG_Class{Session}{Idx+4}(:,k));
                end
            end

            for N = 1:nShuffle
                if rem(N,100) == 0
                    disp(['Currently at ' num2str(N) 'th shuffling of GNG classifier']);
                end
                GNGShuffle{1} = randsample(sum(TrialNum(1:2)),sum(TrialNum(1:2)),'false');
                GNGShuffle{3} = randsample(sum(TrialNum(3:4)),sum(TrialNum(3:4)),'false');
                GNGShuffle{5} = randsample(sum(TrialNum(5:6)),sum(TrialNum(5:6)),'false');
                for Idx = [1 2 3]
                    Trials1 = GNGShuffle{2*Idx-1}(1:TrialNum(2*Idx-1)); Trials2 = GNGShuffle{2*Idx-1}(TrialNum(2*Idx-1)+1:end);
                    for ii = 1:size(CatDecodeResult,1)
                        ShuffleGNGCR{Session}(ii,Idx,N) = (numel(find(CatGNGDecodeResult{Idx}(ii,Trials1) == 2*Idx-1)) + numel(find(CatGNGDecodeResult{Idx}(ii,Trials2) == 2*Idx))) / (TrialNum(2*Idx-1) + TrialNum(2*Idx));
                    end
                end
            end

            Data = MultiTrialScore{S_Idx,Type};
            CatAllInconDecodeResult = [];
            CatAVInconDecodeResult = [];
            CatVAInconDecodeResult = [];
            CatInconVDecodeResult = []; CatInconADecodeResult = [];
            for Idx = 1:4
                InconTrialData{Idx} = Data{Idx+4}(:,PC,:);
                InconTrialNum(Idx) = size(InconTrialData{Idx},3);
                for k = 1:size(InconTrialData{Idx},3)
                    for i = 1:size(InconTrialData{Idx},1)
                        D = [];
                        for j = 1:4
                            D(j) = norm(InconTrialData{Idx}(i,:,k)-AvgIncon{j}(i,:));
                            Inc_DtoUniAvg{S_Idx,Idx}(i,j,k) = norm(InconTrialData{Idx}(i,:,k)-AvgSingle{j}(i,:))/RMS;
                        end
                        [~, InconAllClass{Session}{Idx}(i,k)] = min(D);
                    end
                    CatAllInconDecodeResult = cat(2,CatAllInconDecodeResult,InconAllClass{Session}{Idx}(:,k));
                    for i = 1:size(InconTrialData{Idx},1)
                        D = [];
                        if Idx < 3
                            for j = [1 2]
                                D(j) = norm(InconTrialData{Idx}(i,:,k)-AvgIncon{j}(i,:));
                            end
                        else
                            D(1:2) = nan;
                            for j = [3 4]
                                D(j) = norm(InconTrialData{Idx}(i,:,k)-AvgIncon{j}(i,:));
                            end
                        end
                        [~, InconClass{Session}{Idx}(i,k)] = min(D);
                    end
                    if Idx < 3
                        CatAVInconDecodeResult = cat(2,CatAVInconDecodeResult,InconClass{Session}{Idx}(:,k));
                    else
                        CatVAInconDecodeResult = cat(2,CatVAInconDecodeResult,InconClass{Session}{Idx}(:,k));
                    end
                    for i = 1:size(InconTrialData{Idx},1)
                        D = [];
                        for j = [1 2]
                            D(j) = norm(InconTrialData{Idx}(i,:,k)-AvgSingle{j}(i,:));
                        end
                        [~, InconVClass{Session}{Idx}(i,k)] = min(D);

                        D = []; D(1:2) = nan;
                        for j = [3 4]
                            D(j) = norm(InconTrialData{Idx}(i,:,k)-AvgSingle{j}(i,:));
                        end
                        [~, InconAClass{Session}{Idx}(i,k)] = min(D);
                    end
                    CatInconVDecodeResult = cat(2,CatInconVDecodeResult,InconVClass{Session}{Idx}(:,k));
                    CatInconADecodeResult = cat(2,CatInconADecodeResult,InconAClass{Session}{Idx}(:,k));
                end
            end

            for N = 1:nShuffle
                if rem(N,100) == 0
                    disp(['Currently at ' num2str(N) 'th shuffling of Incon trial']);
                end
                for Idx = 1:4
                    AllShuffle = randsample(sum(InconTrialNum(1:4)),sum(InconTrialNum(1:4)),'false');
                    if Idx == 1
                        Trials = AllShuffle(1:InconTrialNum(Idx));
                    else
                        Trials = AllShuffle(sum(InconTrialNum(1:Idx-1))+1:sum(InconTrialNum(1:Idx)));
                    end
                    for ii = 1:size(CatAVInconDecodeResult,1)
                        AllInconShuffleCR{Session}(ii,Idx,N) = numel(find(CatAllInconDecodeResult(ii,Trials) == Idx)) / InconTrialNum(Idx);
                    end

                    if Idx < 3
                        Shuffle = randsample(sum(InconTrialNum(1:2)),sum(InconTrialNum(1:2)),'false');
                        if Idx == 1
                            Trials = Shuffle(1:InconTrialNum(Idx));
                        else
                            Trials = Shuffle(sum(InconTrialNum(1:Idx-1))+1:sum(InconTrialNum(1:Idx)));
                        end
                        for ii = 1:size(CatAVInconDecodeResult,1)
                            InconShuffleCR{Session}(ii,Idx,N) = numel(find(CatAVInconDecodeResult(ii,Trials) == Idx)) / InconTrialNum(Idx);
                        end
                    else
                        Shuffle = randsample(sum(InconTrialNum(3:4)),sum(InconTrialNum(3:4)),'false');
                        if Idx == 3
                            Trials = Shuffle(1:InconTrialNum(Idx));
                        else
                            Trials = Shuffle(sum(InconTrialNum(3:Idx-1))+1:sum(InconTrialNum(3:Idx)));
                        end
                        for ii = 1:size(CatVAInconDecodeResult,1)
                            InconShuffleCR{Session}(ii,Idx,N) = numel(find(CatVAInconDecodeResult(ii,Trials) == Idx)) / InconTrialNum(Idx);
                        end
                    end
                end
            end
        end
        clearvars PC TrialNum
    end

    switch Type
        case 1; mkdir('Classifier by Stim subspace'); cd('Classifier by Stim subspace');
        case 2; mkdir('Classifier by Deci subspace'); cd('Classifier by Deci subspace');
        case 3; mkdir('Classifier by Indep subspace'); cd('Classifier by Indep subspace');
        case 4; mkdir('Classifier by Interact subspace'); cd('Classifier by Interact subspace');
        otherwise; error('Check type');
    end

    save('Classifier Result.mat');

    %%
    clear variables; close all; clc;

    load('Classifier Result.mat');

    for i = 1:3
        SessionIdx{i} = [];
    end
    SessionIdx{1} = 1:size(Name,1); %mkdir('01 Merged'); mkdir('02 Stationary');  mkdir('03 Moving');
    for i = 1:size(Name,1)
        if strcmp(Name{i,2},'s') == 1; SessionIdx{2} = cat(1,SessionIdx{2},i); end
        if strcmp(Name{i,2},'m') == 1; SessionIdx{3} = cat(1,SessionIdx{3},i); end
    end

    Case = 1;

    TargetIdx = SessionIdx{Case};

    PostTime1 = 0.5;
    PostTime2 = 1;

    DistBin = (Pre+PostTime1)*ImgHz+1:(Pre+PostTime2)*ImgHz;

    clearvars Data DecodeCR CR S_CR AvgCR

    for i = 1:numel(TargetIdx)
        Data = SingleClass{TargetIdx(i)};
        for j = 1:numel(Data)
            for ii = 1:size(Data{j},1)
                CR{j}(ii,i) = numel(find(Data{j}(ii,:) == j))./size(Data{j},2).*100;
            end
            AvgCR{j} = nanmean(CR{j}(DistBin,:),1);
        end
    end

    clearvars S_CI95;
    for i = 1:numel(TargetIdx)
        Data = ShuffleCR{TargetIdx(i)}*100;
        for jj = 1:size(Data,2)
            for ii = 1:size(Data,1)
                S_Avg = nanmean(Data(ii,jj,:),3);
                N = size(Data,3);
                S_Sem = nanstd(Data(ii,jj,:),0,3)./sqrt(N);
                CI95 = tinv([0.025 0.975], N-1);
                yCI95 = bsxfun(@times,S_Sem,CI95(:));
                S_CI95(ii,jj,i) = S_Avg + yCI95(2);
            end
        end
    end

    xAxis = -Pre:1/ImgHz:Post;

    Color = {'r',[243 152 0]./255,'b',[0 152 243]./255};
    Style = {'-','-','-','-'};
    Order = [3 4 1 2];

    MergedBarFig = figure('Position',[50 50 125 110]);
    hold on
    for i = 1:numel(AvgCR)
        Idx = Order(i);

        S_Avg = nanmean(AvgCR{Idx}(SessionIdx{2}));
        S_Sem = nanstd(AvgCR{Idx}(SessionIdx{2}))./sqrt(numel(SessionIdx{2}));
        M_Avg = nanmean(AvgCR{Idx}(SessionIdx{3}));
        M_Sem = nanstd(AvgCR{Idx}(SessionIdx{3}))./sqrt(numel(SessionIdx{3}));

        bar(2*i-1,S_Avg,0.6,'FaceColor','none','lineWidth',1,'EdgeColor',Color{Idx});
        errorbar(2*i-1,S_Avg,S_Sem,S_Sem,'CapSize',3,'color',Color{Idx},'lineWidth',0.75);

        bar(2*i-0.2,M_Avg,0.6,'FaceColor','none','lineWidth',1,'EdgeColor',Color{Idx},'lineStyle',':');
        errorbar(2*i-0.2,M_Avg,M_Sem,M_Sem,'CapSize',3,'color',Color{Idx},'lineWidth',0.75);
    end

    xticks([1 1.8 3 3.8 5 5.8 7 7.8]);
    xticklabels({'S','M','S','M','S','M','S','M'})
    xlim([0.4 8.4]);
    ylim([0 100]);
    set(gca,'TickDir','out','FontName','Arial','FontSize',6);
    xlabel('Type','FontName','Arial','Fontsize',6);
    ylabel('Mean classifier accuracy (%)','FontName','Arial','Fontsize',6);

    saveas(MergedBarFig,'Fig 2l, Stationary Moving Single Classification - Bar.svg');

    SinglePlotFig = figure('Position',[50 50 80 110]);;
    hold on

    for j = 1:size(S_CI95,2)
        plot(xAxis,nanmean(S_CI95(:,SessionIdx{2},:),3),'color',Color{j},'lineWidth',0.5);
    end

    for j = 1:size(S_CI95,2)
        plot(xAxis,nanmean(S_CI95(:,SessionIdx{3},:),3),'color',Color{j},'lineWidth',0.5,'lineStyle',':');
    end

    for j = 1:numel(CR)
        Avg = nanmean(CR{j}(:,SessionIdx{2}),2);
        Sem = nanstd(CR{j}(:,SessionIdx{2}),0,2)./sqrt(numel(SessionIdx{2}));
        errorshade(xAxis,Avg',Avg'+Sem',Avg'-Sem',Color{j});
        alpha(0.2);
        %plot(xAxis,Avg,'color',Color{j},'lineStyle',Style{j},'lineWidth',0.25);
    end

    for j = 1:numel(CR)
        Avg = nanmean(CR{j}(:,SessionIdx{3}),2);
        Sem = nanstd(CR{j}(:,SessionIdx{3}),0,2)./sqrt(numel(SessionIdx{3}));
        errorshade(xAxis,Avg',Avg'+Sem',Avg'-Sem',Color{j});
        alpha(0.2);
        %plot(xAxis,Avg,'color',Color{j},'lineStyle',Style{j},'lineWidth',0.25);
    end

    for j = 1:numel(CR)
        Avg = nanmean(CR{j}(:,SessionIdx{2}),2);
        plot(xAxis,Avg,'color',Color{j},'lineStyle',Style{j},'lineWidth',1);
    end

    for j = 1:numel(CR)
        Avg = nanmean(CR{j}(:,SessionIdx{3}),2);
        plot(xAxis,Avg,'color',Color{j},'lineStyle',Style{j},'lineWidth',1,'lineStyle',':');
    end

    line([0 0],[0 100],'color',[0.75 0.75 0.75],'lineWidth',0.5);
    line([1 1],[0 100],'lineStyle',':','color',[0.75 0.75 0.75],'lineWidth',0.5);
    xlim([-1 1.5]);
    ylim([0 100]);
    set(gca,'TickDir','out','FontName','Arial','FontSize',6);
    xlabel('Time(s)','FontName','Arial','Fontsize',6);
    ylabel('Classifier accuracy (%)','FontName','Arial','Fontsize',6);

    saveas(SinglePlotFig,'Fig 2k, Single Classification - PSTH.svg');

    %
    % StimType{1,1} = 'Vg'; StimType{2,1} = 'Vng'; StimType{3,1} = 'Ag'; StimType{4,1} = 'Ang';
    % LocoType{1,1} = 'S'; LocoType{2,1} = 'M';
    %
    % N = 0; Data = 0;
    % for i = 1:numel(AvgCR)
    %     for ii = 1:numel(SessionIdx{2})
    %         N = N+1;
    %         Data(N,1) = AvgCR{i}(SessionIdx{2}(ii));
    %         Type1{N,1} = StimType{i};
    %         Type2{N,1} = LocoType{1};
    %     end
    %     for ii = 1:numel(SessionIdx{3})
    %         N = N+1;
    %         Data(N,1) = AvgCR{i}(SessionIdx{3}(ii));
    %         Type1{N,1} = StimType{i};
    %         Type2{N,1} = LocoType{2};
    %     end
    % end
    %
    % [p,origtbl,stats] = anovan(Data,{Type1 Type2},"Model","interaction", ...
    %     "Varnames",["Stim","Loco"]);
    %
    % [results,~,~,gnames] = multcompare(stats,"Dimension",[1 2],'CriticalValueType','lsd');
    % %[results,~,~,gnames] = multcompare(stats,"Dimension",[2],'CriticalValueType','hsd');
    %
    % tbl = array2table(results(:,[1 2 6]),"VariableNames", ...
    %     ["Group A","Group B","P-value"]);
    % tbl.("Group A")=gnames(tbl.("Group A"));
    % tbl.("Group B")=gnames(tbl.("Group B"))

    %%
    % Plotting GNG classifier result

    xAxis = -Pre:1/ImgHz:Post;

    Color = {'r','b','m'};
    Style = {'-','-','-'};

    clearvars Data DecodeCR CR S_CR AvgCR

    GNGClassPlotFig = figure('Position',[50 50 80 110]);
    hold on

    TargetIdx = SessionIdx{2};
    for i = 1:numel(TargetIdx)
        Data = GNG_Class{TargetIdx(i)};
        for j = 1:3
            for ii = 1:size(Data{j},1)
                CR{j}(ii,i) = (numel(find(Data{2*j-1}(ii,:) == 2*j-1))+numel(find(Data{2*j}(ii,:) == j*2)))./(size(Data{2*j-1},2)+size(Data{2*j},2)).*100;
            end
            AvgCR{j} = nanmean(CR{j}(DistBin,:),1);
        end
    end

    clearvars S_CI95;
    for i = 1:numel(TargetIdx)
        Data = ShuffleGNGCR{TargetIdx(i)}*100;
        for jj = 1:size(Data,2)
            for ii = 1:size(Data,1)
                S_Avg = nanmean(Data(ii,jj,:),3);
                N = size(Data,3);
                S_Sem = nanstd(Data(ii,jj,:),0,3)./sqrt(N);
                CI95 = tinv([0.025 0.975], N-1);
                yCI95 = bsxfun(@times,S_Sem,CI95(:));
                S_CI95(ii,jj,i) = S_Avg + yCI95(2);
            end
        end
    end

    for j = [1 3]
        plot(xAxis,nanmean(S_CI95(:,j,:),3),'color',Color{j},'lineWidth',0.5);
    end

    for j = [1 3]
        Avg = nanmean(CR{j},2);
        Sem = nanstd(CR{j},0,2)./sqrt(size(CR{ j},2));
        errorshade(xAxis,Avg',Avg'+Sem',Avg'-Sem',Color{j});
        alpha(0.2);
        %plot(xAxis,Avg,'color',Color{j},'lineStyle',Style{j},'lineWidth',0.25);
    end
    for j = [1 3]
        Avg = nanmean(CR{j},2);
        plot(xAxis,Avg,'color',Color{j},'lineStyle',Style{j},'lineWidth',1);
    end

    line([0 0],[0 100],'color',[0.75 0.75 0.75],'lineWidth',0.5);
    line([1 1],[0 100],'lineStyle',':','color',[0.75 0.75 0.75],'lineWidth',0.5);
    xlim([-1 1.5]);
    ylim([0 100]);
    set(gca,'TickDir','out','FontName','Arial','FontSize',6);
    xlabel('Time(s)','FontName','Arial','Fontsize',6);
    ylabel('Classifier accuracy (%)','FontName','Arial','Fontsize',6);

    saveas(GNGClassPlotFig,'Fig S3c, Stationary GNG Classification - PSTH.svg');

    clearvars Data DecodeCR CR S_CR AvgCR

    GNGClassPlotFig = figure('Position',[50 50 80 110]);
    hold on

    TargetIdx = SessionIdx{3};
    for i = 1:numel(TargetIdx)
        Data = GNG_Class{TargetIdx(i)};
        for j = 1:3
            for ii = 1:size(Data{j},1)
                CR{j}(ii,i) = (numel(find(Data{2*j-1}(ii,:) == 2*j-1))+numel(find(Data{2*j}(ii,:) == j*2)))./(size(Data{2*j-1},2)+size(Data{2*j},2)).*100;
            end
            AvgCR{j} = nanmean(CR{j}(DistBin,:),1);
        end
    end

    clearvars S_CI95;
    for i = 1:numel(TargetIdx)
        Data = ShuffleGNGCR{TargetIdx(i)}*100;
        for jj = 1:size(Data,2)
            for ii = 1:size(Data,1)
                S_Avg = nanmean(Data(ii,jj,:),3);
                N = size(Data,3);
                S_Sem = nanstd(Data(ii,jj,:),0,3)./sqrt(N);
                CI95 = tinv([0.025 0.975], N-1);
                yCI95 = bsxfun(@times,S_Sem,CI95(:));
                S_CI95(ii,jj,i) = S_Avg + yCI95(2);
            end
        end
    end

    for j = [1 3]
        plot(xAxis,nanmean(S_CI95(:,j,:),3),'color',Color{j},'lineWidth',0.5);
    end

    for j = [1 3]
        Avg = nanmean(CR{j},2);
        Sem = nanstd(CR{j},0,2)./sqrt(size(CR{ j},2));
        errorshade(xAxis,Avg',Avg'+Sem',Avg'-Sem',Color{j});
        alpha(0.2);
        %plot(xAxis,Avg,'color',Color{j},'lineStyle',Style{j},'lineWidth',0.25);
    end
    for j = [1 3]
        Avg = nanmean(CR{j},2);
        plot(xAxis,Avg,'color',Color{j},'lineStyle',Style{j},'lineWidth',1);
    end

    line([0 0],[0 100],'color',[0.75 0.75 0.75],'lineWidth',0.5);
    line([1 1],[0 100],'lineStyle',':','color',[0.75 0.75 0.75],'lineWidth',0.5);
    xlim([-1 1.5]);
    ylim([0 100]);
    set(gca,'TickDir','out','FontName','Arial','FontSize',6);
    xlabel('Time(s)','FontName','Arial','Fontsize',6);
    ylabel('Classifier accuracy (%)','FontName','Arial','Fontsize',6);

    saveas(GNGClassPlotFig,'Fig S3c, Moving GNG Classification - PSTH.svg');

    GNGBarFig = figure('Position',[50 50 80 110]);
    hold on

    clearvars CR AvgCR

    TargetIdx = SessionIdx{2};
    for i = 1:numel(TargetIdx)
        Data = GNG_Class{TargetIdx(i)};
        for j = 1:3
            for ii = 1:size(Data{j},1)
                CR{j}(ii,i) = (numel(find(Data{2*j-1}(ii,:) == 2*j-1))+numel(find(Data{2*j}(ii,:) == j*2)))./(size(Data{2*j-1},2)+size(Data{2*j},2)).*100;
            end
            AvgCR{j} = nanmean(CR{j}(DistBin,:),1);
        end
    end

    xTick = [3.5 0 1];

    for i = [1 3]
        %Idx = Order(i)
        Avg = nanmean(AvgCR{i});
        Sem = nanstd(AvgCR{i})./sqrt(numel(AvgCR{i}));
        bar(xTick(i),Avg,0.8,'FaceColor','none','lineWidth',1,'EdgeColor',Color{i});
        errorbar(xTick(i),Avg,Sem,Sem,'CapSize',3,'color',Color{i},'lineWidth',0.75);
    end

    TargetIdx = SessionIdx{3};
    for i = 1:numel(TargetIdx)
        Data = GNG_Class{TargetIdx(i)};
        for j = 1:3
            for ii = 1:size(Data{j},1)
                CR{j}(ii,i) = (numel(find(Data{2*j-1}(ii,:) == 2*j-1))+numel(find(Data{2*j}(ii,:) == j*2)))./(size(Data{2*j-1},2)+size(Data{2*j},2)).*100;
            end
            AvgCR{j} = nanmean(CR{j}(DistBin,:),1);
        end
    end

    xTick = [4.5 0 2];

    for i = [1 3]
        %Idx = Order(i)
        Avg = nanmean(AvgCR{i});
        Sem = nanstd(AvgCR{i})./sqrt(numel(AvgCR{i}));
        bar(xTick(i),Avg,0.8,'FaceColor','none','lineWidth',1,'EdgeColor',Color{i});
        errorbar(xTick(i),Avg,Sem,Sem,'CapSize',3,'color',Color{i},'lineWidth',0.75);
    end

    xticks([1 2 3.5 4.5]);
    xticklabels({'S','M','S','M'})
    xlim([0.3 5.2]);
    ylim([0 100]);
    set(gca,'TickDir','out','FontName','Arial','FontSize',6);
    xlabel('Type','FontName','Arial','Fontsize',6);
    ylabel('Mean classifier accuracy (%)','FontName','Arial','Fontsize',6);

    saveas(GNGBarFig,'Fig S3d, GNG Classification - Bar.svg');

    %%
    %plot Incon classifier results

    Data = Uni_DtoUniAvg;
    for i = 1:size(Data,1)
        for j = 1:size(Data,2)
            for ii = 1:size(Data{i,j},1)
                D = []; Label = []; clearvars LDA K L;
                for jj = 1:size(Data{i,j},2)
                    UniStd{j}(ii,jj,i) = std(squeeze(Data{i,j}(ii,jj,:)));
                    D = cat(1,D,squeeze(Data{i,j}(ii,jj,:)));
                    clearvars Temp;
                    if jj == j
                        Temp(1:size(Data{i,j},3),1) = {'Target'};
                        Label = cat(1,Label,Temp);
                    else
                        Temp(1:size(Data{i,j},3),1) = {'NonTarget'};
                        Label = cat(1,Label,Temp);
                    end
                end

                LDA = fitcdiscr(D,Label);

                K = LDA.Coeffs(1,2).Const;
                L = LDA.Coeffs(1,2).Linear;

                Thres(ii,j,i) = abs(K/L);
                % for iii = 1:20
                %     Prc_UniD{iii,1}(ii,j,i) = prctile(squeeze(Data{i,j}(ii,j,:)),5*iii);
                % end
            end
        end
    end
    %
    Data = Inc_DtoUniAvg;
    V_Target = [2 2 1 1];
    for i = 1:size(Data,1)
        for j = 1:size(Data,2)
            for ii = 1:size(Data{i,j},1)
                for k = 1:size(Data{i,j},3)
                    if Data{i,j}(ii,V_Target(j),k) <= Thres(ii,V_Target(j),i);
                        V_DomSort{i,j}(ii,k) = 1;
                    else
                        V_DomSort{i,j}(ii,k) = 0;
                    end
                end
            end
        end
    end
    %
    Data = Inc_DtoUniAvg;
    CorrectAns = [0 1 1 0];
    for i = 1:size(Data,1)
        for S = 1:100
            clearvars Shuffle Idx1 Idx2 Temp
            for j = [1 3]
                Temp = cat(3,Data{i,j},Data{i,j+1});
                Idx1 = randsample(size(Temp,3),size(Data{i,j},3),'false');
                Idx2 = setxor(1:size(Temp,3),Idx1);
                Shuffle{j} = Temp(:,:,Idx1);
                Shuffle{j+1} = Temp(:,:,Idx2);
            end
            for j = 1:4
                for k = 1:size(Shuffle{j},3)
                    for iii = 1:size(Shuffle{j},1)
                        if Shuffle{j}(iii,V_Target(j),k) <= Thres(iii,V_Target(j),i);
                            S_V_DomSort{i,j}(iii,k) = 1;
                        else
                            S_V_DomSort{i,j}(iii,k) = 0;
                        end
                    end
                end
            end
            Case = 0;
            for j = [1 3]
                Case = Case + 1;
                for iii = 1:size(Shuffle{j},1)
                    S_CR{Case}(iii,i,S) = (numel(find(squeeze(S_V_DomSort{i,j}(iii,:)) == CorrectAns(j))) + numel(find(squeeze(S_V_DomSort{i,j+1}(iii,:)) == CorrectAns(j+1))))/(size(S_V_DomSort{i,j},2) + size(S_V_DomSort{i,j+1},2)).* 100;
                end
            end
        end
    end
    %
    Data = V_DomSort;
    CorrectAns = [0 1 1 0];
    for i = 1:size(Data,1);
        Case = 0;
        for j = [1 3];
            Case = Case + 1;
            for ii = 1:size(Data{i,j},1)
                CR{Case}(ii,i) = (numel(find(squeeze(Data{i,j}(ii,:)) == CorrectAns(j))) + numel(find(squeeze(Data{i,j+1}(ii,:)) == CorrectAns(j+1))))/(size(Data{i,j},2)+size(Data{i,j+1},2)).* 100;
            end
        end
    end

    clearvars Color Style;

    Color = [0 0 0; 0 0 0;]./255;
    Style = {':','-'};

    Order = [3.5 1];

    %AnovaTest = [];

    clearvars Avg Sem
    Data = CR;
    for i = 1:2
        Avg(i) = nanmean(nanmean(Data{i}(DistBin,SessionIdx{2}),1)); Sem(i) = nanstd(nanmean(Data{i}(DistBin,SessionIdx{2}),1))./sqrt(numel(SessionIdx{2}));
        %AnovaTest = cat(1,AnovaTest,nanmean(Data{i}(DistBin,SessionIdx{2}),1));
    end

    Fig = figure('Position', [50 50 80 100]);
    hold on;
    for i = 1:2
        bar(Order(i),Avg(i),'FaceColor','none','EdgeColor',Color(i,:),'lineWidth',1,'lineStyle',Style{i});
        errorbar(Order(i),Avg(i),Sem(i),Sem(i),'CapSize',3,'lineWidth',0.75,'color',Color(i,:));
    end

    Order = [4.5 2];

    clearvars Avg Sem
    Data = CR;
    for i = 1:2
        Avg(i) = nanmean(nanmean(Data{i}(DistBin,SessionIdx{3}),1)); Sem(i) = nanstd(nanmean(Data{i}(DistBin,SessionIdx{3}),1))./sqrt(numel(SessionIdx{3}));
        %AnovaTest = cat(1,AnovaTest,nanmean(Data{i}(DistBin,SessionIdx{3}),1));
    end

    for i = 1:2
        bar(Order(i),Avg(i),'FaceColor','none','EdgeColor',Color(i,:),'lineWidth',1,'lineStyle',Style{i});
        errorbar(Order(i),Avg(i),Sem(i),Sem(i),'CapSize',3,'lineWidth',0.75,'color',Color(i,:));
    end

    xticks([1 2 3.5 4.5]);
    xlim([0.3 5.2]);
    xticklabels({'S','M','S','M'});

    set(gca,'TickDir','out','FontName','Arial','FontSize',6);
    ylabel('Mean classification accuracy (%)','FontName','Arial','Fontsize',6);
    xlabel('Type','FontName','Arial','Fontsize',6);

    saveas(Fig,'Fig 3f, Incon Classifier Correct rate - Bar.svg');

    clearvars Avg Sem S_CI95;

    Avg1 = nanmean(CR{2},2); Sem1 = nanstd(CR{2},0,2)./sqrt(size(CR{2},2));
    Avg2 = nanmean(CR{1},2); Sem2 = nanstd(CR{1},0,2)./sqrt(size(CR{1},2));

    xAxis = -Pre:1/ImgHz:Post;

    for j = 1:2
        for jj = 1:size(S_CR{j},2)
            for ii = 1:size(S_CR{j},1)
                Avg = nanmean(S_CR{j}(ii,jj,:),3);
                N = size(S_CR{j},3);
                Sem = nanstd(S_CR{j}(ii,jj,:),0,3)./sqrt(size(S_CR{j},3));
                CI95 = tinv([0.025 0.975], N-1);
                yCI95 = bsxfun(@times,Sem,CI95(:));
                S_CI95{j}(ii,jj) = Avg + yCI95(2);
            end
        end
    end

    xAxis = -1:1/ImgHz:3;

    Fig = figure('Position', [50 50 80 100]);
    hold on;

    errorshade(xAxis,Avg1',Avg1'+Sem1',Avg1'-Sem1',[0 0 0]);
    alpha(0.1);
    errorshade(xAxis,Avg2',Avg2'+Sem2',Avg2'-Sem2',[0 0 0]);
    alpha(0.1);

    plot(xAxis,Avg1,'lineWidth',1,'color','k');
    plot(xAxis,Avg2,'lineWidth',1,'color','k','lineStyle',':');

    plot(xAxis,nanmean(S_CI95{1},2),'lineWidth',0.5,'color','k','lineStyle',':');
    plot(xAxis,nanmean(S_CI95{2},2),'lineWidth',0.5,'color','k');

    line([0 0],[0 100],'lineWidth',0.5,'color',[190 190 190]./255);
    line([1 1],[0 100],'lineWidth',0.5,'color',[190 190 190]./255,'lineStyle',':');

    xlim([-1 1.5]);
    ylim([0 100]);
    set(gca,'TickDir','out','FontName','Arial','FontSize',6);
    xlabel('Time (s)','FontName','Arial','FontSize',6);
    ylabel('Classifier accuracy (%)','FontName','Arial','FontSize',6);

    saveas(Fig,'Fig 3f, Incon Classifier Correct rate - PSTH.svg');

    %%
    %Plot example Incon classification data

    ExIdx = 3;
    OneSecPoint = (Pre+1)*ImgHz+1;

    close all;
    clearvars Color
    Color = [255 0 0; 243 152 0; 0 0 255; 0 152 243]./255;
    %GrayColor = GrayBlackColor(20);
    HistoXBin = 0:1:26;

    Fig = figure('Position',[50 50 110 100]);
    hold on
    for j = 1:4
        histogram(squeeze(Uni_DtoUniAvg{ExIdx,1}(OneSecPoint,j,:)),HistoXBin,'FaceColor',Color(j,:),'lineWidth',0.25);
    end

    line([Thres(OneSecPoint,1,ExIdx) Thres(OneSecPoint,1,ExIdx)],[0 25],'color','r','lineWidth',1,'lineStyle',':');
    xlim([0 20]);
    ylim([0 25]);
    set(gca,'TickDir','out','FontName','Arial','FontSize',6);
    xlabel('Distance to VG trajectory','FontName','Arial','FontSize',6);
    ylabel('#Trial','FontName','Arial','FontSize',6);

    saveas(Fig,'Fig 3e, Distance example - to Vgo.svg');

    Fig = figure('Position',[50 50 110 100]);
    hold on
    for j = 1:1
        histogram(squeeze(Inc_DtoUniAvg{ExIdx,3}(OneSecPoint,j,:)),HistoXBin,'FaceColor',[0 0 0],'lineWidth',0.5);
    end
    for j = 1:1
        histogram(squeeze(Inc_DtoUniAvg{ExIdx,4}(OneSecPoint,j,:)),HistoXBin,'FaceColor',[160 160 160]./255,'lineWidth',0.25);
    end

    line([Thres(OneSecPoint,1,ExIdx) Thres(OneSecPoint,1,ExIdx)],[0 25],'lineWidth',1,'color','r','lineStyle',':');

    xlim([0 20]);
    ylim([0 25]);
    set(gca,'TickDir','out','FontName','Arial','FontSize',6);
    xlabel('Distance to VG trajectory','FontName','Arial','FontSize',6);
    ylabel('#Trial','FontName','Arial','FontSize',6);


    saveas(Fig,'Fig 3e, Distance example - INCON.svg');

    cd ../

    % clearvars Data CR AvgShuffleCR AvgCR
    % Color = {[0.6 0.6 0.6],[0 0 0],[0 0 0],[0.6 0.6 0.6]};
    % Style = {':',':','-','-'};
    %
    % for i = 1:numel(TargetIdx)
    %     Data = InconAllClass{TargetIdx(i)};
    %     if Case == 1
    %         BunZa = zeros(41,1); BunMo = zeros(41,1);
    %     end
    %
    %     for j = 1:numel(Data)
    %         for ii = 1:size(Data{j},1)
    %             InconAllCR{j}(ii,i) = numel(find(Data{j}(ii,:) == j))./size(Data{j},2).*100;
    %             if Case == 1;
    %                 BunZa(ii,1) = BunZa(ii,1) + numel(find(Data{j}(ii,:) == j)); BunMo(ii) = BunMo(ii) + size(Data{j},2);
    %             end
    %         end
    %         AvgCR{j} = nanmean(InconAllCR{j}(DistBin,:),1);
    %         AvgPreCR{j} = nanmean(InconAllCR{j}(1:Pre*ImgHz,:),1);
    %     end
    %     if Case == 1
    %         MergedInconCR(:,i) = BunZa./BunMo.*100;
    %     end
    % end
    % if Case == 1
    %     MergedInconCR_Avg{Type} = nanmean(InconAllCR{j}(DistBin,:),1);
    % end
    %
    % for j = 1:4
    %     S_CR{j} = [];
    % end
    % for i = 1:numel(TargetIdx)
    %     Data = nanmean(AllInconShuffleCR{TargetIdx(i)}.*100,3);
    %     for j = 1:4
    %         S_CR{j} = cat(2,S_CR{j},Data(:,j));
    %     end
    % end
    %
    % Order = [4 1 3 2];
    %
    % InconBarFig = figure('Position',[50 50 80 110]);
    % hold on
    %
    % for i = 1:numel(AvgCR)
    %     Idx = Order(i);
    %     Avg = nanmean(AvgCR{Idx});
    %     Sem = nanstd(AvgCR{Idx})./sqrt(numel(AvgCR{Idx}));
    %     bar(i,Avg,0.6,'FaceColor','none','lineWidth',1,'EdgeColor',Color{Idx},'lineStyle',Style{Idx});
    %     errorbar(i,Avg,Sem,Sem,'CapSize',3,'color',Color{Idx},'lineWidth',0.75);
    % end
    %
    % xticks([1.5 3.5]);
    % xticklabels({'A dom','V dom'})
    % xlim([0.4 4.6]);
    % ylim([0 100]);
    % set(gca,'TickDir','out','FontName','Arial','FontSize',6);
    % xlabel('Type','FontName','Arial','Fontsize',6);
    % ylabel('Mean classifier accuracy (%)','FontName','Arial','Fontsize',6);
    %
    % saveas(InconBarFig,'Classifier 04-2 InconAll_CR_Bar.svg');
    % signrank(AvgCR{1},AvgCR{2})
    % signrank(AvgCR{1},AvgCR{3})
    % signrank(AvgCR{1},AvgCR{4})
    % signrank(AvgCR{2},AvgCR{3})
    % signrank(AvgCR{2},AvgCR{4})
    % signrank(AvgCR{3},AvgCR{4})

    clearvars -except Type
end
