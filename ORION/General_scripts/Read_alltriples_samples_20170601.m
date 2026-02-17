filenames = {'25546ON.csv','25546POST.csv','25546PRE.csv','25979ON.csv','25979POST.csv','25979PRE.csv','26531ON.csv','26531POST.csv','26531PRE.csv','27960ON.csv','27960POST.csv','27960PRE.csv','33466ON.csv','33466POST.csv','33466PRE.csv','33680ON.csv','33680POST.csv','33680PRE.csv','36844ON.csv','36844POST.csv','36844PRE.csv','37648ON.csv','37648POST.csv','37648PRE.csv'};
samplenames ={'sample25546ON','sample25546POST','sample25546PRE','sample25979ON','sample25979POST','sample25979PRE','sample26531ON','sample26531POST','sample26531PRE','sample27960ON','sample27960POST','sample27960PRE','sample33466ON','sample33466POST','sample33466PRE','sample33680ON','sample33680POST','sample33680PRE','sample36844ON','sample36844POST','sample36844PRE','sample37648ON','sample37648POST','sample37648PRE'};
 
allsample = table;

for index = 1:length(filenames)
    myName = filenames{index};
    temp1 = readtable(myName);
    eval(strcat(samplenames{index},'=temp1;'));
    temp1.labels = repmat(samplenames(index),length(temp1.HOECHST1),1);
    if isempty(allsample)
        allsample = temp1;
    else
        allsample = vertcat(allsample,temp1);
    end
end
