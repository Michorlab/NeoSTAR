%%  Calculation entropy for pooled frames
% Jerry lin 2016/12/05


%% Initilization

poolsize = 4;

total_col = 17;
total_row = 13;
cell_limit = 1000;

%cell samplearray;
flag=1;
for i =1:total_col*total_row;
        data1 = cleardata(cleardata.frame == i,:);
        cellno = size(data1,1);
        
        if(cellno>cell_limit)
             samplearray{flag}=data1.EMGM20;
             flag = flag+1;
        end
end

entropy2 = NaN(100,4);

for j=1:100
    sample1 = datasample(samplearray,4);
    sample1=sample1';
    sample1 = cell2mat(sample1);
    
    cellno = size(sample1,1);
    entropy2(j,1) = cellno;
    entropy2(j,2) = 0;
    entropy2(j,3) = 0;
    entropy2(j,4) = 0;
    
    ne1 = nentropy(sample1,'shannon');
    sample2 = datasample(cleardata,cellno);
    ne2 = nentropy(sample2.EMGM20,'shannon');
    entropy2(j,2) = ne1;
    entropy2(j,3) = ne2;
    entropy2(j,4) = ne1/ne2;
end


            