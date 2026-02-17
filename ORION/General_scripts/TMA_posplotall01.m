

data1 = TMA{6};
data2 = data1{:,6:20};
labels = chlabels(6:20);

figure;
for i = 1:3;
    for j = 1:5;
        chs = (i-1)*5+j;
        subplot(3,5,chs);
        scatter(data1.X,data1.Y,log(data2(:,chs)),log(data2(:,chs)),'fill');colormap(jet);
        
        title(labels(chs));
    end
end


      


