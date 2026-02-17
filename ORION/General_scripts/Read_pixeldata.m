%% Read pixel data from ImageJ
%  Jerry Lin 2018/01/26


%% Initialization


myDIR = slideDIR;
myName = slideName;

totalframe = input('Please input freame nubmer:');
mytablename = input('Please input table name:','s'); %eg. 'Vimentin_555_all';

allpixel = cell(length(myName),totalframe);
alltable = table;

%% Read through different cycles

for i=1:length(myName)
    alldata = zeros(totalframe,8);
    
    for j=1:totalframe
        myFile1 = strcat(myDIR{i},'\target-',num2str(j),'.csv');
        myFile2 = strcat(myDIR{i},'\whole-',num2str(j),'.csv');
        %disp(myFile1);
        alldata(j,1) = j;
        if exist(myFile1,'file')
            display(strcat('Processing:',myName{i},'_frame:',num2str(j)));
            
            test1 = readtable(myFile1);
            test2 = readtable(myFile2);
            data1 = test1.Value;
            allpixel{i,j}=data1;
            data2 = test2.Value;
            alldata(j,2) = min(data1);
            alldata(j,3) = max(data1);
            alldata(j,4) = prctile(data1,1);
            alldata(j,5) = prctile(data1,99);
            alldata(j,6) = skewness(data1);
            alldata(j,7) = snr(data1);
            alldata(j,8) = snr(data2);
        end
    end
    
    table1 = array2table(alldata,'VariableNames',{'frame','min','max','pr1','pr99','skewness','snr','snr_total'});
    table1.Cycle = repmat(myName(i),totalframe,1);
    table1 = table1(table1.max>0,:);
    
    if isempty(alltable)
        alltable = table1;
    else
        alltable = vertcat(alltable,table1);
    end
    
        
    clear alldata data1 data2 test1 test2 table1;
end

eval(strcat(mytablename,' = alltable;'));
eval(strcat(mytablename,'_pixel=allpixel;'));


