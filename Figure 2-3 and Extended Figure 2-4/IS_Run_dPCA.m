% This code conducts dPCA on individual data
% You should run this first before testing other dPCA relevant codes.
% Run this code in the folder that has example data

clear variables; close all; clc;

Mat = FindMatFiles();
%clearvars -except EventF

DS_Factor = 2; %Temporal downsampling of data (20Hz to 10Hz)
Pre = 1; %Pre-stimulus time epoch for dPCA
Post = 3; %Post-stimulus time epoch for dCPA

KMTime = 1.5; % Time epoch to dissociate start-moving or keep-moving trial;
KMThreshold = 3; % Speed threshold for keep-moving trial

SmoothFactor = 3; % Size of bins for temporal gaussian smoothing
UseNormalization = 0;

Lick = [1 2; 2 1; 1 2; 2 1;]; %1 = trial with lick (hit or false alarm; 2 = trial with nolick (correct rejection or miss)

for MatIdx = 1:size(Mat)
    
    load(Mat{MatIdx});
    ImgHz = ImgHz/DS_Factor;
    
    nCell = size(EventF{1,1},2); 
    
    Check = 1;
    for i = [1 2 4]
        for j = 1:4
            if isempty(EventF{i,j}) == 1 | size(EventF{i,j},3) < 3;
                Check = 0;
            end
        end
    end
    
    if Check == 1
        TypeIdx = 0;
        for i = 1:4
            for j = 1:4
                FData = []; DS_FData = []; CutEventF = [];
                TypeIdx = TypeIdx+1;
                if isempty(EventF{i,j}) == 0
                    if DS_Factor == 1
                        DS_FData = EventF{i,j};
                    else
                        for kk = 1:size(EventF{i,j},3); for jj = 1:nCell; for ii = 1:size(EventF{i,j},1)/DS_Factor;
                                    DS_FData(ii,jj,kk) = nanmean(EventF{i,j}(DS_Factor*(ii-1)+1:DS_Factor*ii,jj,kk),1);
                        end; end; end;
                    end

                    if SmoothFactor > 0
                        DS_FData = smoothdata(DS_FData,1,'gaussian',SmoothFactor);
                    end
                    
                    for kk = 1:size(DS_FData,3)
                        CutEventF(:,:,kk) = DS_FData((PSTHPre-Pre)*ImgHz:(PSTHPre+Post)*ImgHz,:,kk);
                    end
                    
                    InitialTrialNum(TypeIdx) = size(CutEventF,3);
                    ArrangedF{TypeIdx} = CutEventF;
                    ArrangeSpeed{TypeIdx} = EventSpeed{i,j};
                    for k = 1:size(ArrangedF{TypeIdx},3)
                        Avg = nanmean(ArrangedF{TypeIdx}(1:Pre*ImgHz,:,k),1);
                        Std = nanstd(ArrangedF{TypeIdx}(1:Pre*ImgHz,:,k),0,1);
                        for ii = 1:size(ArrangedF{TypeIdx},1)
                            NormArrangedF{TypeIdx}(ii,:,k) = (ArrangedF{TypeIdx}(ii,:,k)-Avg);
                        end
                    end
                else
                    NormArrangedF{TypeIdx} = [];
                    InitialTrialNum(TypeIdx) = 0;
                    ArrangedF{TypeIdx} = [];
                end
            end
        end
        
        N = nCell;
        S = 4;
        T = (Pre+Post)*ImgHz;
        D = 2;
        E = max(InitialTrialNum);
        
        if UseNormalization == 1
            ArrangedF = NormArrangedF;
        end
        
        clearvars StimF StimSpeed
        
        for i = 1:4
            for j = [1 3]
                if j == 1
                    StimF{i,1} = cat(3,ArrangedF{4*(i-1)+j},ArrangedF{4*(i-1)+(j+1)});
                    StimSpeed{i,1} = cat(1,EventSpeed{i,j},EventSpeed{i,j+1});
                elseif j == 3
                    StimF{i,2} = cat(3,ArrangedF{4*(i-1)+j},ArrangedF{4*(i-1)+(j+1)});
                    StimSpeed{i,2} = cat(1,EventSpeed{i,j},EventSpeed{i,j+1});
                end
            end
        end
        for i = 1:4
            InconF{i} = ArrangedF{12+i};
            InconSpeed{i} = EventSpeed{4,i};
        end
        % firingRates: N x S x D x T x maxTrialNum
        for n = 1:nCell
            for s = 1:S
                for t = 1:size(ArrangedF{s*2-1},1)
                    if isempty(ArrangedF{s*2-1}) == 0
                        for k = 1:size(ArrangedF{s*2-1},3)
                            FiringRates(n,s,Lick(s,1),t,k) = ArrangedF{s*2-1}(t,n,k);
                            TrialNum(n,s,Lick(s,1)) = size(ArrangedF{s*2-1},3);
                        end
                    else
                        TrialNum(n,s,Lick(s,1)) = 0;
                    end
                end
                for t = 1:size(ArrangedF{s*2},1)
                    if isempty(ArrangedF{s*2}) == 0
                        for k = 1:size(ArrangedF{s*2},3)
                            FiringRates(n,s,Lick(s,2),t,k) = ArrangedF{s*2}(t,n,k);
                            TrialNum(n,s,Lick(s,2)) = size(ArrangedF{s*2},3);
                        end
                    else
                        TrialNum(n,s,Lick(s,2)) = 0;
                    end
                end
            end
        end
        
        for n = 1:nCell; for s = 1:S; for d = 1:D;
                    FiringRates(n,s,d,:,TrialNum(n,s,d)+1:E) = nan;
        end; end; end;
        
        % firingRates: N x S x D x T x maxTrialNum
        for n = 1:nCell
            for s = 1:S
                for t = 1:size(ArrangedF{(s+4)*2-1},1)
                    if isempty(ArrangedF{(s+4)*2-1}) == 0
                        for k = 1:size(ArrangedF{(s+4)*2-1},3)
                            MultiFiringRates(n,s,Lick(s,1),t,k) = ArrangedF{(s+4)*2-1}(t,n,k);
                            MultiTrialNum(n,s,Lick(s,1)) = size(ArrangedF{(s+4)*2-1},3);
                        end
                    elseif isempty(ArrangedF{(s+4)*2-1}) == 1
                        MultiTrialNum(n,s,Lick(s,1)) = 0;
                    end
                end
                for t = 1:size(ArrangedF{(s+4)*2},1)
                    if isempty(ArrangedF{(s+4)*2}) == 0
                        for k = 1:size(ArrangedF{(s+4)*2},3)
                            MultiFiringRates(n,s,Lick(s,2),t,k) = ArrangedF{(s+4)*2}(t,n,k);
                            MultiTrialNum(n,s,Lick(s,2)) = size(ArrangedF{(s+4)*2},3);
                        end
                    elseif isempty(ArrangedF{(s+4)*2}) == 1
                        MultiTrialNum(n,s,Lick(s,2)) = 0;
                    end
                end
            end
        end
        
        for n = 1:nCell
            for s = 1:S
                for d = 1:D
                    MultiFiringRates(n,s,d,:,MultiTrialNum(n,s,d)+1:E) = nan;
                end
            end
        end
        
        FiringRatesAverage = nanmean(FiringRates, 5);
        MultiFiringRatesAverage = nanmean(MultiFiringRates, 5);
        
        if strcmp(State,'s') == 1
            for i = 1:4
                for j = 1:2
                    S_FiringRates{i,j} = []; S_FiringRatesAvg{i,j} = [];
                    for n = 1:nCell
                        for t = 1:size(StimF{i,j},1);
                            for k = 1:size(StimF{i,j},3)
                                S_FiringRates{i,j}(n,1,1,t,k) = StimF{i,j}(t,n,k);
                            end
                        end
                    end
                    S_FiringRatesAvg{i,j} = nanmean(S_FiringRates{i,j},5);
                end
            end

            for jj = 1:4
                if isempty(InconF{jj})
                    InconTrialNumber{MatIdx,1}(3,jj) = 0;
                else
                    InconTrialNumber{MatIdx,1}(3,jj) = size(InconF{jj},3);
                end

                S_InconF = InconF{jj};
                for n = 1:nCell
                    for t = 1:size(S_InconF,1)
                        for k = 1:size(S_InconF,3)
                            S_InconFiringRates{jj}(n,1,1,t,k) = S_InconF(t,n,k);
                        end
                    end
                end
                S_InconFiringRatesAvg{jj} = nanmean(S_InconFiringRates{jj},5);
            end
        elseif strcmp(State,'m') == 1
            for i = 1:4
                for j = 1:2
                    KMIdx = []; SMIdx = [];
                    KM_FiringRates{i,j} = []; KM_FiringRatesAvg{i,j} = [];
                    SM_FiringRates{i,j} = []; SM_FiringRatesAvg{i,j} = [];
                    for k = 1:size(StimSpeed{i,j},1)
                        if mean(StimSpeed{i,j}(k,1:KMTime*EventHz)) > KMThreshold
                            KMIdx = cat(1,KMIdx,k);
                        else
                            SMIdx = cat(1,SMIdx,k);
                        end
                    end
                    KM_StimF = StimF{i,j}(:,:,KMIdx);
                    SM_StimF = StimF{i,j}(:,:,SMIdx);
                    for n = 1:nCell
                        for t = 1:size(KM_StimF,1)
                            for k = 1:size(KM_StimF,3)
                                KM_FiringRates{i,j}(n,1,1,t,k) = KM_StimF(t,n,k);
                            end
                            for k = 1:size(SM_StimF,3)
                                SM_FiringRates{i,j}(n,1,1,t,k) = SM_StimF(t,n,k);
                            end
                        end
                    end
                    KM_FiringRatesAvg{i,j} = nanmean(KM_FiringRates{i,j},5);
                    SM_FiringRatesAvg{i,j} = nanmean(SM_FiringRates{i,j},5);
                end
            end
            
            for jj = 1:4
                KMIdx = []; SMIdx = [];
                KM_InconFiringRates{jj} = []; KM_InconFiringRatesAvg{jj} = [];
                SM_InconFiringRates{jj} = []; SM_InconFiringRatesAvg{jj} = [];
                for k = 1:size(InconSpeed{jj},1)
                    if mean(InconSpeed{jj}(k,1:KMTime*EventHz)) > KMThreshold
                        KMIdx = cat(1,KMIdx,k);
                    else
                        SMIdx = cat(1,SMIdx,k);
                    end
                end
                
                InconTrialNumber{MatIdx,1}(1,jj) = numel(KMIdx);
                InconTrialNumber{MatIdx,1}(2,jj) = numel(SMIdx);
                
                KM_InconF = InconF{jj}(:,:,KMIdx);
                SM_InconF = InconF{jj}(:,:,SMIdx);
                for n = 1:nCell
                    for t = 1:size(KM_InconF,1)
                        for k = 1:size(KM_InconF,3)
                            KM_InconFiringRates{jj}(n,1,1,t,k) = KM_InconF(t,n,k);
                        end
                        for k = 1:size(SM_InconF,3)
                            SM_InconFiringRates{jj}(n,1,1,t,k) = SM_InconF(t,n,k);
                        end
                    end
                end
                KM_InconFiringRatesAvg{jj} = nanmean(KM_InconFiringRates{jj},5);
                SM_InconFiringRatesAvg{jj} = nanmean(SM_InconFiringRates{jj},5);
            end
        end
        
        ifSimultaneousRecording = true;
        
        CombinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
        MargNames = {'Stimulus', 'Decision', 'Condition-independent', 'S-D Interaction'};
        MargColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
        
        Time = -Pre:1/ImgHz:Post;
        TimeEvents = 0;

        %%%%%%%%%%%% Step 1: PCA of the dataset
        X = FiringRatesAverage(:,:);
        [W,~,~] = svd(X, 'econ');
        if nCell < 20
            W = W(:,1:nCell);
        else
            W = W(:,1:20);
        end
        %minimal plotting
        dpca_plot(FiringRatesAverage, W, W, @dpca_plot_default);
        
        % computing explained variance
        explVar = dpca_explainedVariance(FiringRatesAverage, W, W, ...
            'combinedParams', CombinedParams);
        
        % a bit more informative plotting
        dpca_plot(FiringRatesAverage, W, W, @dpca_plot_default, ...
            'explainedVar', explVar, ...
            'time', Time,                        ...
            'timeEvents', TimeEvents,               ...
            'marginalizationNames', MargNames, ...
            'marginalizationColours', MargColours);
        
        %%%%%%%%%%%% Step 2: PCA in each marginalization separately
        dpca_perMarginalization(FiringRatesAverage, @dpca_plot_default, ...
            'combinedParams', CombinedParams);
        
        %%%%%%%%%%%% Step 3: dPCA without regularization and ignoring noise covariance
        
        % This is the core function.
        % W is the decoder, V is the encoder (ordered by explained variance),
        % whichMarg is an array that tells you which component comes from which
        % marginalization

        tic
        if nCell > 20
            [W,V,whichMarg] = dpca(FiringRatesAverage, 20, ...
                'combinedParams', CombinedParams);
        else
            [W,V,whichMarg] = dpca(FiringRatesAverage, nCell, ...
                'combinedParams', CombinedParams);
        end
        toc

        explVar = dpca_explainedVariance(FiringRatesAverage, W, V, ...
            'combinedParams', CombinedParams);

        dpca_plot(FiringRatesAverage, W, V, @dpca_plot_default, ...
            'explainedVar', explVar, ...
            'marginalizationNames', MargNames, ...
            'marginalizationColours', MargColours, ...
            'whichMarg', whichMarg,                 ...
            'time', Time,                        ...
            'timeEvents', TimeEvents,               ...
            'timeMarginalization', 3, ...
            'legendSubplot', 16);

        %%%%%%%%%%%% Step 4: dPCA with regularization

        % This function takes some minutes to run. It will save the computations
        % in a .mat file with a given name. Once computed, you can simply load
        % lambdas out of this file:
        %   load('tmp_optimalLambdas.mat', 'optimalLambda')

        % Please note that this now includes noise covariance matrix Cnoise which
        % tends to provide substantial regularization by itself (even with lambda set
        % to zero).


        optimalLambda = dpca_optimizeLambda(FiringRatesAverage, FiringRates, TrialNum, ...
            'combinedParams', CombinedParams, ...
            'simultaneous', ifSimultaneousRecording, ...
            'numRep', 10, ...  % increase this number to ~10 for better accuracy
            'filename', 'tmp_optimalLambdas.mat');

        Cnoise = dpca_getNoiseCovariance(FiringRatesAverage, ...
            FiringRates, TrialNum, 'simultaneous', ifSimultaneousRecording);
        if nCell > 20
            [W,V,whichMarg] = dpca(FiringRatesAverage, 20, ...
                'combinedParams', CombinedParams, ...
                'lambda', optimalLambda, ...
                'Cnoise', Cnoise);
        else
            [W,V,whichMarg] = dpca(FiringRatesAverage, nCell, ...
                'combinedParams', CombinedParams, ...
                'lambda', optimalLambda, ...
                'Cnoise', Cnoise);
        end

        if nCell < 20
            W = W(:,1:nCell);
            V = V(:,1:nCell);
        end

        explVar = dpca_explainedVariance(FiringRatesAverage, W, V, ...
            'combinedParams', CombinedParams);

        dpca_plot(FiringRatesAverage, W, V, @dpca_plot_default, ...
            'explainedVar', explVar, ...
            'marginalizationNames', MargNames, ...
            'marginalizationColours', MargColours, ...
            'whichMarg', whichMarg,                 ...
            'time', Time,                        ...
            'timeEvents', TimeEvents,               ...
            'timeMarginalization', 3,           ...
            'legendSubplot', 16);

        XFull = FiringRatesAverage;
        X = XFull(:,:)';
        DataDim = size(XFull);
        Z = X * W;

        XFull_Test = MultiFiringRatesAverage;
        X_Test = XFull_Test(:,:)';
        TestDataDim = size(XFull_Test);
        TestZ = X_Test * W;
        
        NumOfStimuli = 4;
        
        UniColors = {[1 0 0],[1 0.5 0],[0 0 1],[0 0.5 1]};
        MultiColors = {[0.5 0.5 0.5],[0 0 0];[0 0 0],[0.5 0.5 0.5]};
        MultiStyles = {':',':';'-','-'};
        
        Title = split(Mat{MatIdx},'_');
        Title = strjoin(Title);
        
        xAxis = -Pre:1/ImgHz:Post;
        
        if MatIdx == 1; mkdir('dPCA Result'); end; cd('dPCA Result');
        
        for MargIdx = 1:4
            
            Idx = find(whichMarg  == MargIdx);

            for i = 1:numel(Idx)
                
                componentsToPlot = Idx(i);
                Data = reshape(Z(:,componentsToPlot)', [length(componentsToPlot) DataDim(2:end)]);
                
                for f = 1:NumOfStimuli
                    AvgScore{MatIdx,MargIdx}{2*f-1}(:,i) = squeeze(Data(1,f,1,:));
                    AvgScore{MatIdx,MargIdx}{2*f}(:,i) = squeeze(Data(1,f,2,:));
                    
                    % firingRates: N x S x D x T x maxTrialNum
                    for k = 1:TrialNum(1,f,1)
                        X_SingleTrial = FiringRates(:,f,1,:,k);
                        X_Trial = X_SingleTrial(:,:)';
                        Z_Trial = X_Trial * W;
                        SingleDim = size(X_SingleTrial);
                        Z_TrialFull = reshape(Z_Trial(:,componentsToPlot)', [length(componentsToPlot) SingleDim(2:end)]);
                        UniTrialScore{MatIdx,MargIdx}{2*f-1}(:,i,k) = squeeze(Z_TrialFull(1,1,1,:));
                    end
                    for k = 1:TrialNum(1,f,2)
                        X_SingleTrial = FiringRates(:,f,2,:,k);
                        X_Trial = X_SingleTrial(:,:)';
                        Z_Trial = X_Trial * W;
                        SingleDim = size(X_SingleTrial);
                        Z_TrialFull = reshape(Z_Trial(:,componentsToPlot)', [length(componentsToPlot) SingleDim(2:end)]);
                        UniTrialScore{MatIdx,MargIdx}{2*f}(:,i,k) = squeeze(Z_TrialFull(1,1,1,:));
                    end
                end
                
            end
            
            for  i = 1:numel(Idx)
                
                componentsToPlot = Idx(i);
                TestData = reshape(TestZ(:,componentsToPlot)', [length(componentsToPlot) TestDataDim(2:end)]);
                
                for f = 1:NumOfStimuli
                    AvgScore{MatIdx,MargIdx}{2*(f+4)-1}(:,i) = squeeze(TestData(1,f,1,:));
                    AvgScore{MatIdx,MargIdx}{2*(f+4)}(:,i) = squeeze(TestData(1,f,2,:));
                    
                    for k = 1:MultiTrialNum(1,f,1)
                        X_MultiTrial = MultiFiringRates(:,f,1,:,k);
                        X_Trial = X_MultiTrial(:,:)';
                        Z_Trial = X_Trial * W;
                        MultiDim = size(X_MultiTrial);
                        Z_TrialFull = reshape(Z_Trial(:,componentsToPlot)', [length(componentsToPlot) MultiDim(2:end)]);
                        MultiTrialScore{MatIdx,MargIdx}{2*f-1}(:,i,k) = squeeze(Z_TrialFull(1,1,1,:));
                    end
                    for k = 1:MultiTrialNum(1,f,2)
                        X_MultiTrial = MultiFiringRates(:,f,2,:,k);
                        X_Trial = X_MultiTrial(:,:)';
                        Z_Trial = X_Trial * W;
                        MultiDim = size(X_MultiTrial);
                        Z_TrialFull = reshape(Z_Trial(:,componentsToPlot)', [length(componentsToPlot) MultiDim(2:end)]);
                        MultiTrialScore{MatIdx,MargIdx}{2*f}(:,i,k) = squeeze(Z_TrialFull(1,1,1,:));
                    end
                end
            end
            
            if strcmp(State,'s') == 1
                for i = 1:4
                    for j = 1:2
                        FullTestData = S_FiringRatesAvg{i,j};
                        TestData = FullTestData(:,:)';
                        TestDataResult = TestData * W;
                        for ii = 1:numel(Idx)
                            MvmtSortedScore{MatIdx,MargIdx}{i,j}{3}(:,ii) = TestDataResult(:,Idx(ii));
                        end
                        
                        InputData = S_FiringRates{i,j};
                        for k = 1:size(InputData,5)
                            FullTestData = InputData(:,:,:,:,k);
                            TestData = FullTestData(:,:)';
                            TestDataResult = TestData * W;
                            for ii = 1:numel(Idx)
                                MvmtSortedTrialScore{MatIdx,MargIdx}{i,j}{3}(:,ii,k) = TestDataResult(:,Idx(ii));
                            end
                        end
                    end
                end
                
                for jj = 1:4
                    if isempty(S_InconFiringRatesAvg{jj}) == 0
                        FullTestData = S_InconFiringRatesAvg{jj};
                        TestData = FullTestData(:,:)';
                        TestDataResult = TestData * W;
                        for ii = 1:numel(Idx)
                            InconMvmtSortedScore{MatIdx,MargIdx}{1,jj}{3}(:,ii) = TestDataResult(:,Idx(ii));
                        end
                    else
                        InconMvmtSortedScore{MatIdx,MargIdx}{1,jj}{3} = [];
                    end
                end
            elseif strcmp(State,'m') == 1
                for i = 1:4
                    for j = 1:2
                        FullTestData = SM_FiringRatesAvg{i,j};
                        TestData = FullTestData(:,:)';
                        TestDataResult = TestData * W;
                        for ii = 1:numel(Idx)
                            MvmtSortedScore{MatIdx,MargIdx}{i,j}{2}(:,ii) = TestDataResult(:,Idx(ii));
                        end
                        
                        InputData = SM_FiringRates{i,j};
                        for k = 1:size(InputData,5)
                            FullTestData = InputData(:,:,:,:,k);
                            TestData = FullTestData(:,:)';
                            TestDataResult = TestData * W;
                            for ii = 1:numel(Idx)
                                MvmtSortedTrialScore{MatIdx,MargIdx}{i,j}{2}(:,ii,k) = TestDataResult(:,Idx(ii));
                            end
                        end
                        
                        FullTestData = KM_FiringRatesAvg{i,j};
                        TestData = FullTestData(:,:)';
                        TestDataResult = TestData * W;
                        for ii = 1:numel(Idx)
                            MvmtSortedScore{MatIdx,MargIdx}{i,j}{1}(:,ii) = TestDataResult(:,Idx(ii));
                        end
                        
                        InputData = KM_FiringRates{i,j};
                        for k = 1:size(InputData,5)
                            FullTestData = InputData(:,:,:,:,k);
                            TestData = FullTestData(:,:)';
                            TestDataResult = TestData * W;
                            for ii = 1:numel(Idx)
                                MvmtSortedTrialScore{MatIdx,MargIdx}{i,j}{1}(:,ii,k) = TestDataResult(:,Idx(ii));
                            end
                        end
                    end
                end
                for jj = 1:4
                    if isempty(SM_InconFiringRatesAvg{jj}) == 0
                        FullTestData = SM_InconFiringRatesAvg{jj};
                        TestData = FullTestData(:,:)';
                        TestDataResult = TestData * W;
                        for ii = 1:numel(Idx)
                            InconMvmtSortedScore{MatIdx,MargIdx}{1,jj}{2}(:,ii) = TestDataResult(:,Idx(ii));
                        end
                    else
                        InconMvmtSortedScore{MatIdx,MargIdx}{1,jj}{2} = [];
                    end
                    
                    if isempty(KM_InconFiringRatesAvg{jj}) == 0
                        FullTestData = KM_InconFiringRatesAvg{jj};
                        TestData = FullTestData(:,:)';
                        TestDataResult = TestData * W;
                        for ii = 1:numel(Idx)
                            InconMvmtSortedScore{MatIdx,MargIdx}{1,jj}{1}(:,ii) = TestDataResult(:,Idx(ii));
                        end
                    else
                        InconMvmtSortedScore{MatIdx,MargIdx}{1,jj}{1} = [];
                    end
                end
            end
        end
        
        explVar.whichMarg = whichMarg;
        ComponentInfo{MatIdx,1} = explVar;
        Name{MatIdx,1} = Mat{MatIdx};
        Name{MatIdx,2} = State;
        Conv.V{MatIdx} = V;
        Conv.W{MatIdx} = W;
        NumCell(MatIdx,1) = nCell;
        
        cd ../
    end
    
    clearvars -except Mat MatIdx DS_Factor Pre Post UseNormalization Score ComponentInfo ImgHz AvgMultiScore SingleScore Lick SmoothFactor KMTime KMThreshold MvmtSortedScore AvgScore UniTrialScore MultiTrialScore Name InconMvmtSortedScore MvmtSortedTrialScore InconTrialNumber kkk MargNames Conv NumCell
    close all

    delete tmp_optimalLambdas.mat
end

cd('dPCA Result');

%%
save('dPCA Result.mat');


