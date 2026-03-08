load figS3f

%%
CNO1 = conAll(:,:,2);
CNO2 = conAll(:,:,4);

% remove very low values 
CNO1adj = nan(size(CNO1,1),size(CNO1,2));
for iRow = 1:size(CNO1,1)
    for iCol = 1:size(CNO1,2)
        if CNO1(iRow,iCol) > 0.05
            CNO1adj(iRow,iCol) = CNO1(iRow,iCol);
        end
    end
end

CNO2adj = nan(size(CNO2,1),size(CNO2,2));
for iRow = 1:size(CNO2,1)
    for iCol = 1:size(CNO2,2)
        if CNO2(iRow,iCol) > 0.05
            CNO2adj(iRow,iCol) = CNO2(iRow,iCol);
        end
    end
end
        
for iCluster = 1:length(CNO1adj)
    [CNO1sort(iCluster,:),fieldRanks] = sort(CNO1adj(iCluster,:),'descend','missingPlacement','last');
    CNO2sort(iCluster,:) = CNO2adj(iCluster,fieldRanks);
end

CNO1rank = nan(length(CNO1adj),6,5);
CNO2rank = nan(length(CNO2adj),6,5);

for iCluster = 1:length(CNO1adj)
    for compField = 1:5
        fieldIdx = 1 + compField;
        for iField = fieldIdx:6
            CNO1rank(iCluster,iField-1,compField) = calc.diffScore(CNO1sort(iCluster,iField),CNO1sort(iCluster,compField));
        end
    end
end

for iCluster = 1:length(CNO2adj)
    for compField = 1:5
        fieldIdx = 1 + compField;
        for iField = fieldIdx:6
            CNO2rank(iCluster,iField-1,compField) = calc.diffScore(CNO2sort(iCluster,iField),CNO2sort(iCluster,compField));
        end
    end
end

CNO1rankCol = reshape(CNO1rank,[],1);
CNO2rankCol = reshape(CNO2rank,[],1);

x = CNO1rankCol;
y = CNO2rankCol;
xnan = isnan(x);
ynan = isnan(y);
x = x(~xnan & ~ynan);
y = y(~xnan & ~ynan);
xMax = nanmax(x);
yMax = nanmax(y);

figure;
hold on
plot(x,y,'k.','markersize',20, 'Color', '#bdbdbd')
set(gca,'fontsize',14,'fontweight','bold')
axis([0 1 -1 1])


r2 = calc.benFit(x,y,'degree',1,'showLine',1)
[r,p] = corr(x,y)
sum(~isnan(x))

xlabel('Subfield relationships CNO1')
ylabel('Subfield relationships CNO2')

text(0.4,-0.7,sprintf('Control: r = %.2f, p = %.3e',r, p),'color','#969696','fontweight','bold','fontsize',14)

%% hM3
CNO1 = expAll(:,:,2);
CNO2 = expAll(:,:,4);

% remove very low values 
CNO1adj = nan(size(CNO1,1),size(CNO1,2));
for iRow = 1:size(CNO1,1)
    for iCol = 1:size(CNO1,2)
        if CNO1(iRow,iCol) > 0.05
            CNO1adj(iRow,iCol) = CNO1(iRow,iCol);
        end
    end
end

CNO2adj = nan(size(CNO2,1),size(CNO2,2));
for iRow = 1:size(CNO2,1)
    for iCol = 1:size(CNO2,2)
        if CNO2(iRow,iCol) > 0.05
            CNO2adj(iRow,iCol) = CNO2(iRow,iCol);
        end
    end
end
        
for iCluster = 1:length(CNO1adj)
    [CNO1sort(iCluster,:),fieldRanks] = sort(CNO1adj(iCluster,:),'descend','missingPlacement','last');
    CNO2sort(iCluster,:) = CNO2adj(iCluster,fieldRanks);
end

CNO1rank = nan(length(CNO1adj),6,5);
CNO2rank = nan(length(CNO2adj),6,5);

for iCluster = 1:length(CNO1adj)
    for compField = 1:5
        fieldIdx = 1 + compField;
        for iField = fieldIdx:6
            CNO1rank(iCluster,iField-1,compField) = calc.diffScore(CNO1sort(iCluster,iField),CNO1sort(iCluster,compField));
        end
    end
end

for iCluster = 1:length(CNO2adj)
    for compField = 1:5
        fieldIdx = 1 + compField;
        for iField = fieldIdx:6
            CNO2rank(iCluster,iField-1,compField) = calc.diffScore(CNO2sort(iCluster,iField),CNO2sort(iCluster,compField));
        end
    end
end

CNO1rankCol = reshape(CNO1rank,[],1);
CNO2rankCol = reshape(CNO2rank,[],1);

x = CNO1rankCol;
y = CNO2rankCol;
xnan = isnan(x);
ynan = isnan(y);
x = x(~xnan & ~ynan);
y = y(~xnan & ~ynan);
xMax = nanmax(x);
yMax = nanmax(y);

figure;
hold on
plot(x,y,'k.','markersize',20, 'Color', '#525252')
set(gca,'fontsize',14,'fontweight','bold')
axis([0 1 -1 1])

r2 = calc.benFit(x,y,'degree',1,'showLine',1) 
[r,p] = corr(x,y)
sum(~isnan(x))

xlabel('Subfield relationships CNO1')
ylabel('Subfield relationships CNO2')

text(0.4,-0.6,sprintf('hM3: r = %.2f, p = %.3e',r, p),'color','k','fontweight','bold','fontsize',14)