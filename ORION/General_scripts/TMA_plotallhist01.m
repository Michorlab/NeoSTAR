%% Plot all histogram 

figure;

for i=7:24;
    subplot(3,6,i-6);
    histfit(log(alldata{:,i}),100,'kernel');
    title(labels(i));
end
