%%% 20.02.17
%%% This program aim to draw and realign the line profile on 2D images
%%% It needs tif of a single plan with multiple channel
%%% The choosen channel is displayed to draw every profile line needed
%%% The line has to be drawn from left to right !!!
%%% Image will then be rotated to have an horizontal line. A window around the line is defined
%%% Range of adjustable size are used to bin the data to improve realignment. 
%%% Realignent in done based on the selected channel. The profile along the
%%% selected line is then performed and saved in the Result variable.
%%% 
%%% The program give the posibilty to select only a range of data according
%%% to the intensity level of one of the channel, or range of population. The result is then
%%% plotted and fitted with a gaussian curve.

%%% Adjustable parameter : Staining, Window for line profile, width of the
%%% line for realignment, range for intensity selection, Range around the junction to visualize (L197) 
%%% Not yet : /
%%% Save for different range of intensity

clc;
clear;
close all;
%% Charge images
%path='/Users/Xumei/Dropbox/MBI PhD/ExptData/Stretching/Bipul201611/';
path='/Users/nus/Desktop/Tiff-ZO1-Occlu-Ecad/linescan_data';

folders=uipickfiles('FilterSpec',path);
disp([num2str(length(folders)), ' images found']);
Passage=1;
%% Test line profile
resultsfile=fopen(fullfile(folders{1},'Summary_Data.txt'),'w');
fprintf(resultsfile,'Image \tStaining \tValue\n');
Compteur=1;
RoiCell=zeros(1,length(folders));
%%
for NumCell=1:length(folders);
    files=dir(folders{NumCell});
    ind=cellfun(@(s) strfind(s,'.tif'),{files(:).name},'UniformOutput',false);
    filestif=files(~cellfun(@isempty,ind));
    ind2=cellfun(@(s) strfind(s,'.roi'),{files(:).name},'UniformOutput',false);
    roi_file=files(~cellfun(@isempty,ind2));
    filename=strcat(fullfile(folders{NumCell},filestif.name));
    if Passage
        Nchannel=inputdlg('Number of channel ?');
        Nchannel=str2double(Nchannel{1});
    end
    ImCell=uint16(mtifread4D(filename,Nchannel));
    % If first passage
    if Passage
    %Info for each channel
        for leg=1:size(ImCell,4)
            Text=strcat('Protein for channel ',num2str(leg));
            str{leg}=char(inputdlg(Text,'Name the different channel',1));
        end
    end
%    %Channel selection
%     [channel,v] = listdlg('PromptString','Channel for line selection:',...
%                     'SelectionMode','single',...
%                     'ListString',cellstr(str));
%     end
%     %% Line profile until the window with MSG is closed
%     figure('Name','Line Profile');
%     % /!\ Always select line Left to Right
%     h2 = msgbox('Close this window when last profile');
%     %Compteur for number of line draw
%     Count=0;
%     while ishandle(h2),
%         imshow(ImCell(:,:,channel),[min(min(ImCell(:,:,channel))) max(max(ImCell(:,:,channel)))]) 
%         set(gcf, 'Position', get(0, 'Screensize'));
%         if ishandle(h2);
%             h1 = imline();    
%             pos = getPosition(h1);
%             Count=Count+1;
%             LineCoor(Count,1:2)=ceil(pos(1,:));
%             LineCoor(Count,3:4)=ceil(pos(2,:));
%         end
%         clc;
%         disp([num2str(Count), ' lines draw']);
%         close all;
%     end
    LineCoor=zeros(length(roi_file),4);
    for NRoi=1:length(roi_file)
        roifilename=strcat(fullfile(folders{NumCell},roi_file(NRoi).name));
        roitemp=ReadImageJROI(roifilename);
        LineCoor(NRoi,:)=roitemp.vnLinePoints;
    end
    RoiCell(NumCell)=NRoi;
    channel=roitemp.vnPosition(1);
    %% Image rotation for horizontal line
    for j=1:NRoi
        roifilename=strcat(fullfile(folders{NumCell},roi_file(j).name));
        roitemp=ReadImageJROI(roifilename);        
        zselect=roitemp.vnPosition(2);
        % Slope
        slope = (LineCoor(j,2)-LineCoor(j,4))/(LineCoor(j,1)-LineCoor(j,3));

        % Rotation angle
        w = atand(slope);
        RotateImage=imrotate(ImCell(:,:,zselect,:),w,'bilinear','crop');
        RotateImage=squeeze(RotateImage);
        % New line coordinate /!\ Be sure original image is an even-square
        rads=deg2rad(w);
        center=size(ImCell(:,:,1))/2+.5;
        
        y1= (LineCoor(j,2)-center(1))*cos(rads)-(LineCoor(j,1)-center(2))*sin(rads);
        x1=(LineCoor(j,2)-center(1))*sin(rads)+(LineCoor(j,1)-center(2))*cos(rads);
        y1=round(y1+center(1));
        x1=round(x1+center(2));
        
        y2= (LineCoor(j,4)-center(1))*cos(rads)-(LineCoor(j,3)-center(2))*sin(rads);
        x2=(LineCoor(j,4)-center(1))*sin(rads)+(LineCoor(j,3)-center(2))*cos(rads);
        y2=round(y2+center(1));
        x2=round(x2+center(2));   
        
        
        %% Test rotation
%         ImageTest=RotateImage(:,:,channel);
%         ImageBF=ImCell(:,:,zselect,:);
%         figure
%         subplot(1,2,1)
%         imshow(ImageBF(:,:,channel),[min(min(ImCell(:,:,channel))) max(max(ImCell(:,:,channel)))]) 
%         hold on
%         plot(LineCoor(j,1),LineCoor(j,2),'r.','MarkerSize',20)
%         plot(LineCoor(j,3),LineCoor(j,4),'r.','MarkerSize',20)
%         hold off
%         subplot(1,2,2)
%         imshow(ImageTest,[min(min(ImCell(:,:,channel))) max(max(ImCell(:,:,channel)))]) 
%         hold on
%         plot(x1,y1,'r.','MarkerSize',20)
%         plot(x2,y2,'r.','MarkerSize',20) 
%         hold off
        %% Window size Selection
        AreaTemp=20;
        if y1+AreaTemp<size(RotateImage,1) & y1-AreaTemp>0 & x2<size(RotateImage,1) & x1>0 & x1<size(RotateImage,1)% If line selection not to close to border
            Ok=0;
            if x1>x2
                temp=x2;
                x2=x1;
                x1=temp;
            end
            if Passage
                default={'20'};
                Area=20;
                while Ok==0,
                    Imoy=double(RotateImage(:,:,channel));
                    %for display: Imoy en RGB
                    ImoyRGB=zeros([size(Imoy) 3]);
                    Itemp=Imoy/max(Imoy(:));
                    ImoyRGB(:,:,1)=Imoy/max(Imoy(:));
                    ImoyRGB(:,:,2)=Imoy/max(Imoy(:));
                    Itemp(y1-Area:y2+Area,x1:x2)=1;
                    ImoyRGB(:,:,3)=Itemp;
                    hseg=figure('Name','Width selection');
                    figure(hseg)
                    imshow(ImoyRGB),title('Selection')
                    set(gcf, 'Position', get(0, 'Screensize'));
                    text2=['Line width selection ok ? Current ',num2str(Area),'px'];
                    %Test if selection is ok
                    choice = questdlg(text2, ...
                        'Selection', ...
                        'Yes','No','Yes');
                    % Handle response
                    switch choice
                        case 'Yes'
                            Ok=1;
                        case 'No'
                            answer = inputdlg('New window size','Selection',1,default);
                            Area=str2double(answer{:});
                            close all;
                    end
                end
            end
            close all;
            if y1+Area<size(RotateImage,1) & y1-Area>0 & x2<size(RotateImage,1) & x1>0% If line selection not to close to border
                disp(['Cell',num2str(NumCell), ' junction ',num2str(j)]);

        %         %%
                Imcrop=RotateImage(y1-Area:y2+Area,x1:x2,:);
        %         figure;
        %         imshow(Imcrop(:,:,channel),[min(min(Imcrop(:,:,channel))) max(max(Imcrop(:,:,channel)))]) 
        % 

                %% Mean profile along the junction

                Profiletemp=mean(Imcrop,1);
                Profiletemp2=zeros(size(Profiletemp));
                Bkgrd1=sum(Profiletemp(1,1:5,:))/5;
                Bkgrd2=sum(Profiletemp(1,end-4:end,:))/5;
                Bkgrd=squeeze((Bkgrd1+Bkgrd2)/2);
                for f=1:size(Profiletemp,3)
                    Profiletemp2(1,:,f)=Profiletemp(1,:,f)-Bkgrd(f);
                end

                %Realignment based on expected junction size

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ProfileSelection=30;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                DistX1X2=x2-x1;
                if DistX1X2>ProfileSelection*2;  % If Line is longer than the ProfileSelection*2
                    IndexRealign=find(Profiletemp2(1,:,channel)==max(Profiletemp2(1,:,channel),[],2));
                    if IndexRealign==ProfileSelection
                        IndexRealign=ProfileSelection+1;
                    end
                    Profile=Profiletemp2(1,IndexRealign-ProfileSelection:IndexRealign+(ProfileSelection-1),:);
                else
                    IndexRealign=find(Profiletemp2(1,:,channel)==max(Profiletemp2(1,:,channel),[],2));
                    if IndexRealign==ProfileSelection
                        Profile=zeros(1,ProfileSelection*2,size(Profiletemp,3));
                        Profile(1,1:ProfileSelection-IndexRealign+DistX1X2+1,:)=Profiletemp2(1,:,:);        
                    elseif IndexRealign>ProfileSelection               
                        Profile=zeros(1,ProfileSelection*2,size(Profiletemp,3));
 %                       Profile(1,IndexRealign-ProfileSelection:IndexRealign-ProfileSelection+DistX1X2,:)=Profiletemp2(1,:,:);
                        if size(Profile,2)-IndexRealign+ProfileSelection~=size(Profiletemp2,2)-IndexRealign+ProfileSelection+1
                            Diff=(size(Profiletemp2,2)-IndexRealign+ProfileSelection+1)-(size(Profile,2)-IndexRealign+ProfileSelection);
                            if Diff>0;
                            Profile(1,1:end-(IndexRealign-ProfileSelection),:)=Profiletemp2(1,IndexRealign-ProfileSelection+Diff:end,:);  
                            else
                            Profile(1,1:end-(IndexRealign-ProfileSelection-Diff),:)=Profiletemp2(1,IndexRealign-ProfileSelection:end,:);  
    
                            end
                        else
                            Profile(1,1:end-(IndexRealign-ProfileSelection),:)=Profiletemp2(1,IndexRealign-ProfileSelection:end,:);  
                        end
                    else               
                        Profile=zeros(1,ProfileSelection*2,size(Profiletemp,3));
                        Profile(1,ProfileSelection-IndexRealign:ProfileSelection-IndexRealign+DistX1X2,:)=Profiletemp2(1,:,:);        
                    end
                end
                if Passage
                Result=Profile;
                else
                Result=cat(1,Result,Profile(1,1:ProfileSelection*2,:));
                end
                Passage=0;        
            end
        end
        %%
    end
    
end
%% Save results
for ch=1:size(Result,3) 
    NextCell=1;
    NFolders=1;
    for cell=1:size(Result,1)
        if cell==NextCell
            NextCell=NextCell+RoiCell(NFolders);            
            fprintf(resultsfile,[folders{NFolders},'\t']);
            NFolders=NFolders+1;
        else
            fprintf(resultsfile,'\t');
        end
        fprintf(resultsfile,str{ch});
        for res=1:size(Result,2)
            fprintf(resultsfile,['\t',num2str(Result(cell,res,ch))]);
        end
        fprintf(resultsfile,'\n');
    end
end    
%% Figure
figfinal=figure;
ChNumber=Nchannel;
Average=mean(Result,1);
hold on
Name2=('Summary.fig');
for Ncell=1:size(Result,1)
    subplot(1,ChNumber,1)
    hold on
    plot(Result(Ncell,:,1),'--b')
    plot(Average(:,:,1),'-k','LineWidth',2)    
    legend(str{1},'Average')
    hold off
    if ChNumber>1
        subplot(1,ChNumber,2)
        hold on
        plot(Result(Ncell,:,2),'--g') 
        plot(Average(:,:,2),'-k','LineWidth',2)
        legend(str{2},'Average')    
        hold off
    end
    if ChNumber>2
        subplot(1,ChNumber,3)
        hold on
        plot(Result(Ncell,:,3),'--r')   
        plot(Average(:,:,3),'-k','LineWidth',2)
        legend(str{3},'Average')
        hold off
    end
end
saveas(figfinal,fullfile(folders{1},Name2));
% Average
for ch=1:size(Profile,3)
    fprintf(resultsfile,['Average\t',str{ch}]);
    for res=1:size(Profile,2)
        fprintf(resultsfile,['\t',num2str(Average(:,res,ch))]);
    end
    fprintf(resultsfile,'\n');

end

fclose(resultsfile);

%% Homemade Fit
%Gaussian function + offset
%a1*exp(-((x-b1)/c1)^2)+d1
%% 
Average=permute(Average,[2 1 3]);
GaussianOffset = fittype(@(a, b, c, d, x, y) a*exp(-((x-b)/c).^2)+d, 'coefficients', {'a', 'b', 'c', 'd'}, 'independent', {'x'}, 'dependent', {'y'});
x=1:size(Average,1);
h=figure('Name','All data fitted');
MaxAve=max(Average);
for ch=1:size(Average,3)
    [ResultFit,gof] = fit(x', Average(:,:,ch), GaussianOffset, ...
        'StartPoint', [MaxAve(ch), ProfileSelection, 10, 50], ...
        'Lower', [0, ProfileSelection-10, 0, 0], ...
        'Upper', [Inf, ProfileSelection+10, 100, Inf], ...
        'Robust', 'LAR' );
% Figure
    figure(h)
    subplot(1,size(Average,3),ch)
    plot(Average(:,:,ch),'-b')
    hold on
    plot(ResultFit,'-r')
    hold off
    legend(str{ch},'Fit')
    text('units', 'normalized','Position',[0.7 0.86],'Color','k','String',...
    sprintf('Fit results :\n Amp=%g \n width=%g \n offset=%g',ResultFit.a,ResultFit.c*2.355,ResultFit.d))
%     %Width Half Max = 2.355*c
end
set(gcf, 'Position', get(0, 'Screensize'));
saveas(h,fullfile(folders{1},'MeanGauss_Fit.fig'));

%%%%%%%%%%%%

%%%%%%%%%%%%

%%%%%%%%%%%%

%%%%%%%%%%%%

%%%%%%%%%%%%

%%%%%%%%%%%%

%%%%%%%%%%%%

%%%%%%%%%%%%

%% Segmentation based on population pourcentage 
figure(figfinal)
set(gcf, 'Position', get(0, 'Screensize'));

%%
choice2 = questdlg('Data selection ?', ...
                    'Selection', ...
                    'Yes','No','Yes');
            switch choice2
        case 'Yes'
            resultsfileFiltre=fopen(fullfile(folders{1},'Summary_DataSelected.txt'),'w');
            fprintf(resultsfileFiltre,'Junction \tStaining \tValue\n');            
            [ChSelect,v] = listdlg('PromptString','Channel for segmentation ?',...
                'SelectionMode','single',...
                'ListString',cellstr(str));
            choice3 = questdlg('Segmentation by Population or Intensity?', ...
                    'Selection', ...
                    'Population','Intensity','Population');
            switch choice3
                case 'Population'    
                    prompt = {'Min (%):','Max (%):'};
                    dlg_title = 'Population pourcentage';
                    num_lines = 1;
                    defaultans = {'0','100'};
                    range = inputdlg(prompt,dlg_title,num_lines,defaultans);
                    %%
                    BorneMin=round(str2double(range{1})*size(Result,1)/100);
                    BorneMax=round(str2double(range{2})*size(Result,1)/100);
                    if BorneMin==0
                        BorneMin=1;
                    end
                    Resulttemp=permute(Result,[2 1 3]);
                    ResultSelect=Result(:,:,ChSelect);
                    MatrixAverage=max(ResultSelect,[],2);
                    [values,order]=sort(MatrixAverage(:,1));
                    IndexSelec=order(BorneMin:BorneMax);
                case 'Intensity'
                    prompt = {'Min expression level (%):','Max expression level (%):'};
                    dlg_title = 'Expression level';
                    num_lines = 1;
                    defaultans = {'0','100'};
                    range = inputdlg(prompt,dlg_title,num_lines,defaultans);
                    Resulttemp=permute(Result,[2 1 3]);
                    ValMax=max(max(Resulttemp(:,:,ChSelect),[],1));
                    % Find where maxima is in the selected range
                    IndexSelec=find(max(Resulttemp(:,:,ChSelect),[],1)/ValMax*100>=str2double(range{1}) & max(Resulttemp(:,:,ChSelect),[],1)/ValMax*100<=str2double(range{2}));
            end

%%
            ResultFiltre=Resulttemp(:,IndexSelec,:);
            % Result copy
            Compteur=1;
%             for line=1:size(ResultFiltre,2)
%                 fprintf(resultsfileFiltre,num2str(IndexSelec(Compteur)));
%                 Compteur=Compteur+1;
%                 for ch=1:size(ResultFiltre,3)
%                     fprintf(resultsfileFiltre,['\t',str{ch}]);
%                     for res=1:size(Profile,1)
%                         fprintf(resultsfileFiltre,['\t',num2str(ResultFiltre(res,line,ch))]);
%                     end
%                     fprintf(resultsfileFiltre,'\n');
%                 end
%             end
%%          
            ResultFiltre=permute(ResultFiltre,[2 1 3]);
            for ch=1:size(ResultFiltre,3) 
                NextCell=1;
                NFolders=1;
                for cell=1:size(ResultFiltre,1)
                    if cell==NextCell
                        NextCell=NextCell+RoiCell(NFolders);            
                        fprintf(resultsfile,[folders{NFolders},'\t']);
                        NFolders=NFolders+1;
                    else
                        fprintf(resultsfileFiltre,'\t');
                    end
                    fprintf(resultsfileFiltre,str{ch});
                    for res=1:size(Result,2)
                        fprintf(resultsfileFiltre,['\t',num2str(ResultFiltre(cell,res,ch))]);
                    end
                    fprintf(resultsfileFiltre,'\n');
                end
            end 
            %%
            % Final figure
            figfinal=figure;
            ChNumber=size(ResultFiltre,3);
            AverageFiltre=mean(ResultFiltre,1);
            hold on
            Name2=('SummarySelected.fig');
            for Ncell=1:size(ResultFiltre,1)
                subplot(1,ChNumber,1)
                hold on
                plot(ResultFiltre(Ncell,:,1),'--b')
                plot(AverageFiltre(:,:,1),'-k','LineWidth',2)    
                legend(str{1},'Average')
                hold off
                if ChNumber>1
                    subplot(1,ChNumber,2)
                    hold on
                    plot(ResultFiltre(Ncell,:,2),'--g') 
                    plot(AverageFiltre(:,:,2),'-k','LineWidth',2)
                    legend(str{2},'Average')    
                    hold off
                end
                if ChNumber>2
                    subplot(1,ChNumber,3)
                    hold on
                    plot(ResultFiltre(Ncell,:,3),'--r')   
                    plot(AverageFiltre(:,:,3),'-k','LineWidth',2)
                    legend(str{3},'Average')
                    hold off
                end
            end

            saveas(figfinal,fullfile(folders{1},Name2));
            %% Save Average in txt
            AverageFiltre=permute(AverageFiltre,[2 1 3]);
            for ch=1:size(ResultFiltre,3)
                fprintf(resultsfileFiltre,['Average\t',str{ch}]);
                for res=1:size(ResultFiltre,2)
                    fprintf(resultsfile,['\t',num2str(AverageFiltre(res,:,ch))]);
                end
                fprintf(resultsfileFiltre,'\n');
            end
            % Homemade Fit
            GaussianOffset = fittype(@(a, b, c, d, x, y) a*exp(-((x-b)/c).^2)+d, 'coefficients', {'a', 'b', 'c', 'd'}, 'independent', {'x'}, 'dependent', {'y'});
            x=1:size(AverageFiltre,1);
            h=figure('Name','Selected data fitted');
            MaxAve=max(AverageFiltre);
            for ch=1:size(AverageFiltre,3)
                [ResultFit,gof] = fit(x', AverageFiltre(:,:,ch), GaussianOffset, ...
                    'StartPoint', [MaxAve(ch), 30, 10, 50], ...
                    'Lower', [0, 0, 0, 0], ...
                    'Upper', [Inf, 200, 100, Inf], ...
                    'Robust', 'LAR' );
            % Figure
                figure(h)
                subplot(1,size(AverageFiltre,3),ch)
                plot(AverageFiltre(:,:,ch),'-b')
                hold on
                plot(ResultFit,'-r')
                hold off
                legend(str{ch},'Fit')
                text('units', 'normalized','Position',[0.7 0.86],'Color','k','String',...
                sprintf('Fit results :\n Amp=%g \n width=%g \n offset=%g',ResultFit.a,ResultFit.c*2.355,ResultFit.d))
%                 %Width Half Max = 2.355*c
            end
            set(gcf, 'Position', get(0, 'Screensize'));
            saveas(h,fullfile(folders{1},'MeanGauss_FitSelected.fig'));
            fclose(resultsfileFiltre);
                case 'No'
            end
                     
            

