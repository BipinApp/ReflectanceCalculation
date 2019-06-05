clc
clear 
startTime = datetime('now');
% directory scan in matlab to search file  
dirName = 'Z:\ImageDrive\OLI.TIRS\L8\P181\R040\';

%get the directory list
dir_list = dir(dirName);
%checking for the folder 
%first 2 data returned by directory doesn't contains any folder information
% removing the first two rows
dir_list = dir_list(~ismember({dir_list.name},{'.','..'}));
numberOfImages = length(dir_list);

% directory with no LC1 directory 
%preallocation of the vector
%emptyDir = zeros(1,20);
mean_value = zeros(numberOfImages,7);
std_value = zeros(numberOfImages,7);
date_of_year  = transpose(zeros(1,numberOfImages));
date_vector = transpose(zeros(1,numberOfImages));
bandColors={'b','c','g','r','m','[0.6667 0 0]','k'};
bandName={'CA','Blue','Green','Red','NIR','SWIR1','SWIR2'};
bands = 7;
SAA_ROI = zeros(1,numberOfImages);
SZA_ROI = zeros(1,numberOfImages);
VAA_ROI = zeros(1,numberOfImages);
VZA_ROI = zeros(1,numberOfImages);

% store empty acquisition for removal
empty_acquisition = [];

% looping through the directory to search for all the  image aquisiton datas  
for dateAqui= 1:numberOfImages
    %selecting the LC1 folder inside selected aquisition time
    % LC1 is the combined and the lasted data available
    imageLocation = fullfile(dir_list(dateAqui).folder,dir_list(dateAqui).name,'LC1');
    
    % print to terminal to see the processing
    warning('Directory Processing : %s',imageLocation);
    
    if exist(imageLocation,'dir')~= 7
       warning('on')
       warningText = strcat('Folder not Found', imageLocation);
       warning(warningText);
       %emptyDir(dateAqui) = 1;
       empty_acquisition = [empty_acquisition, dateAqui];
       continue;
    end 
    
    %searching the directory for the common string in the file
    file_list = dir(imageLocation);
    
    %reading the common string in the filename
    %filename has constant character length  so this method works here
    % need to change this code to find a common method
    baseName = file_list(3,:).name(1:40);
    
    %clearing the file_list 
    clear file_list
    
    %reading the mtl file using the built in function name MTL_parser_L8
    mtlFileName = fullfile(imageLocation,strcat(baseName,'_','MTL.txt'));
    if exist(mtlFileName,'file') ~=2
        warningText = strcat('MTL file not Found ',mtlFileName);
        warning(warningText);
        empty_acquisition = [empty_acquisition, dateAqui];
        continue
    end
    mtl_list = MTL_parser_L8(mtlFileName);
    
    %reading the acquisition date from the MTL file
    dateOfAqusition = mtl_list.PRODUCT_METADATA.DATE_ACQUIRED;
    % date2doy function is modified for this purpose
    [doy,fraction,dp] = date2doy(datenum(dateOfAqusition));
    date_vector(dateAqui) = dp;
    
    %reading the Band quality assesment file (BQA) file 
    % https://landsat.usgs.gov/collectionqualityband
    bqaFileName = fullfile(imageLocation,strcat(baseName,'_','BQA.TIF'));
    [bqaFile,R] = geotiffread(bqaFileName);
    [bqaRowSize,bqaColSize] = size(bqaFile);
    
    %FILTERING THE IMAGE 
    % LANDSAT 8 COLLECTION 1QA VALUE INTERPRETATION
    % 2720 cirrus confidence low, snow/Ice confidence low, cloud shadow
    % confidece low, cloud No, Radiometric Saturation No, Terrain occlusion
    % No and Fill No  and clear terrain
    % Mask to remove the clouds 
    primaryMask =(bqaFile == 2720);
    
    % correction for the solar zenith and azimuth 
    % Landsat data for the for both solar and view information of the image
    % multibanread -> file consist of two dimensional data
    % 1 -> solar azimuth Angel  and 2 -> solar zenith Angle 
    solarFileName = fullfile(imageLocation,strcat(baseName,'_solar_B05.img'));
    solarAngleInfoAZen = multibandread(solarFileName,[size(bqaFile),2],...
        'int16',0,'bsq','ieee-le');
    
    solarAzimuth = solarAngleInfoAZen(:,:,1)./100;
    solarZenith  = solarAngleInfoAZen(:,:,2)./100; 
    
    % 1 -> sensor azimuth Angel  and 2 -> sensor zenith Angle 
    sensorFileName =fullfile(imageLocation,strcat(baseName,'_sensor_B05.img'));
    sensorAngleInfoAZen = multibandread(sensorFileName,[size(bqaFile),2],...
        'int1 6',0,'bsq','ieee-le');
    
    sensorAzimuth = sensorAngleInfoAZen(:,:,1)./100;
    sensorZenith = sensorAngleInfoAZen(:,:,2)./100;
    
    clear solarAngleInfoAZen solarFileName sensorAngleInfoAZen sensorFileName ...
        
    
    %defining Image ROI
    % Libya4 ROI corner co-ordinated used by CNES
    upperLeftLong = 723825; upperLeftLati = 3171375;
    upperRightLong = 743355; upperRightLati = 3171825;
    lowerLeftLong = 724245; lowerLeftLati = 3149325;
    lowerRightLong = 743805; lowerRightLati = 3149685;
    
    % repeating the co-ordinates to close the ROI
    % chaging the co-ordinates position to get Rectangular ROI
    xCordi = [upperLeftLong upperRightLong lowerRightLong lowerLeftLong  upperLeftLong];
    yCordi = [upperLeftLati upperRightLati lowerRightLati lowerLeftLati  upperLeftLati];
    
    % converting the map co-ordinates to pixel co-ordinates
    [rowPixel, colPixel] = map2pix(R,xCordi,yCordi);
    
    %selection of the region of interest (ROI selection)
    roiSel = poly2mask(rowPixel,colPixel,bqaRowSize,bqaColSize);
    
    % Generating mask based on roiSel and Primasy Mask
    cloudAndCirrusMask = primaryMask .* roiSel;
    solarAzimuthMask = solarAzimuth .* cloudAndCirrusMask;
    solarZenithMask = solarZenith .* cloudAndCirrusMask;
    sensorAzimuthMask = sensorAzimuth.* cloudAndCirrusMask;
    sensorZenithMask = sensorZenith .* cloudAndCirrusMask;
    
    % removing the zero value in the  Mask
    solarAzimuthMask(solarAzimuthMask == 0) = NaN;
    solarZenithMask(solarZenithMask == 0) = NaN;
    sensorAzimuthMask(sensorAzimuthMask == 0) = NaN;
    sensorZenithMask(sensorZenithMask == 0) = NaN;
    
    % calculating the mean value of the angle
    %{
         https://landsat.usgs.gov/what-are-orbit-paths-landsat-satellites
    %}
    SAA_ROI(dateAqui) = nanmean(solarAzimuthMask(:));
    SZA_ROI(dateAqui) = nanmean(solarZenithMask(:));
    VAA_ROI(dateAqui) = nanmean(sensorAzimuthMask(:));
    VZA_ROI(dateAqui) = nanmean(sensorZenithMask(:));
    
    % mtl_parser needs update. 
    % Didn't update because every one is using the same file
    % some information like band information isn't included in the file 
    % some error like spaces in the band in mtl file
    % hardcoding the bands in the code.
    % reflective bands are from 1-7,
    for band = 1:bands
        % reading the specific band for calculation
        bandFileName = fullfile(imageLocation,strcat(baseName,'_B',num2str(band),'.TIF'));
        [bandImage,R] = geotiffread(bandFileName);
        %calculation of TOA reflectance 
        % formula to TOA Reflectance =>  rho = gain*DN + bias
        % Reflectance ranges from 0 to 1, change image type from unsigned 16 bit
        % inter to double
        multiBand = strcat('REFLECTANCE_MULT_BAND_',num2str(band));
        addBand = strcat('REFLECTANCE_ADD_BAND_',num2str(band));
        TOAReflectance = (mtl_list.RADIOMETRIC_RESCALING.(sprintf('%s',multiBand))*double(bandImage)+...
         mtl_list.RADIOMETRIC_RESCALING.(sprintf('%s',addBand)))./cosd(solarZenith);
     
        % converting the roiSel to double to mask the TOA reflectance
        maskedTOAReflectance = TOAReflectance.*double(roiSel);
 
        % masking multiplies by 0 and 1 value
        % converting 0 value to NaN value to exclude in computation of mean and
        % standard deviation
        maskedTOAReflectance(maskedTOAReflectance == 0) = NaN;
        
        % nanmean and nanstd removes NaN value from computation
        mean_value(dateAqui,band) = nanmean(maskedTOAReflectance(:));
        std_value(dateAqui,band) = nanstd(maskedTOAReflectance(:));
        
        % saving the date information for each mean and std-deviation
        if band == 7
            mean_value(dateAqui,band+1) = str2double(dir_list(dateAqui).name);
            mean_value(dateAqui,band+2) = dp;
            std_value(dateAqui,band+1) = str2double(dir_list(dateAqui).name);
            std_value(dateAqui,band+2) = dp;
        end
    end 
    warning('Processing of above directory Complete')
end

% calculating the total time taken to process the data
% takes about one hour to process the data
stopTime = datetime('now');
timeTaken = stopTime - startTime;
warning('Data Processing Complete!!!')

% removing the empty folder directory
mean_value(empty_acquisition,:) = [];
std_value(empty_acquisition,:) = [];
SAA_ROI(empty_acquisition) = [];
SZA_ROI(empty_acquisition) = [];
VAA_ROI(empty_acquisition) = [];
VZA_ROI(empty_acquisition) = [];


% saving the workspace 
save ('mean_value_181_40.mat','mean_value','std_value','empty_acquisition',...
    'bandColors','bandName','SAA_ROI','SZA_ROI','VAA_ROI','VZA_ROI');

%plotting the mean and standard deviation without filtering any data
% plotting temporal uncertainity 
figure
for band = 1:bands
   e = plot(mean_value(:,9),mean_value(:,band));
   e.Marker = 'o';
   e.MarkerSize = 10;
   e.Color = bandColors{band};
   e.LineWidth = 0.75;
   e.MarkerFaceColor = bandColors{band};
   e.MarkerEdgeColor = bandColors{band};
   % removes the connecting line between the points
   e.LineStyle = 'None';
   xlabel('Year')
   ylabel('TOA Reflectance','FontSize',16,'FontWeight','bold')
   title('TOA Reflectance of Libya 4 (Landsat 8) ','FontSize',16,'FontWeight','bold')
   legend(bandName{1:7},'location','northeast')
   hold on
   grid on
   % limit the Y-axis 
   %ylim ([0, 1])
end

%plotting the temporal uncertainity with spatial uncertainiy.
figure 
for band = 1:bands
   e = errorbar(mean_value(:,9),mean_value(:,band),std_value(:,band));
   e.Marker = 'o';
   e.MarkerSize = 10;
   e.Color = bandColors{band};
   e.LineWidth = 0.75;
   e.CapSize = 10;
   e.LineStyle = 'None';
   e.MarkerFaceColor = bandColors{band};
   e.MarkerEdgeColor = bandColors{band};
   xlabel('Year')
   ylabel('TOA Reflectance','FontSize',16,'FontWeight','bold')
   title('TOA Reflectance (Unfiltered Data)','FontSize',16,'FontWeight','bold')
   legend(bandName{1:7},'location','northeast')
   hold on
   grid on
end

%{
All the outlier removal techniques defines a threshold around the data. 
This boundary around the data underestimates the drift in the sensor.
Visual inspection through ENVI to identify cloudy scene is more
appropriate.
Discarding the Technique defined below and defining the cloudy scene
observed in the scence collected.
%}
stdOfMean = std(mean_value(:,1:bands)); 
meanOfMean = mean(mean_value(:,1:bands));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    OUTLIERS REMOVAL Techniques
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The empirical rule tells you what percentage of your data falls within a
%  certain number of standard deviations from the mean:
% 68% of the data falls within 1 standard deviation of the mean.
% 95% of the data falls within 2 standard deviations of the mean.
% 99.7% of the data falls within 3 standard deviations of the mean. 

% PROPERTIES OF NORMAL DISTRIBUTION

% The mean, mode and median are all equal.
% The curve is symmetric at the center (i.e. around the mean, ?).
% Exactly half of the values are to the left of center 
% The total area under the curve is 1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% preallocation of the variables
err.outliers = [];
err.outliersG = [];
err.outliersL = [];

% defining Sigma Value for filtering the data
% used to define the filter and plot titleText
sigma  = 3;

% searching for the outliers
for band = 1:bands
    % greater than outliers 3 \sigma
    err.outliersG = [err.outliersG, find(mean_value(:,band)>...
        (meanOfMean(band) + sigma*stdOfMean(band)))'];
    
    % Less than outliers 3 \sigma
    err.outliersL = [err.outliersL, find(mean_value(:,band)<...
        (meanOfMean(band) - sigma*stdOfMean(band)))'];
end 
% ouliers contains the repeated data so unique is used to remove duplicates
err.outliers = unique([err.outliersG, err.outliersL]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



err.outliers = [60,65,74,83,86,88]; % visual inspection of the Data

% assigning new variable for filtering 
filteredMean = mean_value;
filteredStd = std_value;

%removal of the outliers from the filtered data
filteredMean(err.outliers,:) = [];
filteredStd(err.outliers,:) = [];
SAA_ROI(err.outliers) = [];
SZA_ROI(err.outliers) = [];
VAA_ROI(err.outliers) = [];
VZA_ROI(err.outliers) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        BRDF CORRECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}

% polar plot
figure
sza = polarplot(deg2rad(SZA_ROI),filteredMean(:,1),'DisplayName','Solar Zenith Angle');
sza.LineStyle = 'none';
sza.Marker = '*';
sza.Color = 'RED';
hold on
saa = polarplot(deg2rad(SAA_ROI),filteredMean(:,1),'DisplayName','Solar Azimuth Angle');
saa.LineStyle = 'none';
saa.Marker = 'p';
saa.Color = 'Blue';
hold on
vza= polarplot(deg2rad(VZA_ROI),filteredMean(:,1),'DisplayName','View Zenith Angle');
vza.LineStyle = 'none';
vza.Marker = 'X';
vza.Color = 'RED';
hold on
vaa= polarplot(deg2rad(VAA_ROI),filteredMean(:,1),'DisplayName','View Azimuth Angle');
vaa.LineStyle = 'none';
vaa.Marker = 'h';
vaa.Color = 'Blue';
hold off
lgd = legend('show');
title(lgd, 'Angle')


% converting the polar co-ordinates to cartesian co-ordinates
x_co_solar  = sind(SZA_ROI).*cosd(SAA_ROI);
y_co_solar = sind(SZA_ROI).*sind(SAA_ROI);
z_co_solar = cosd(SZA_ROI); ... Gives height info ... Not included in BRDF

x_co_view = sind(VZA_ROI).*cosd(VAA_ROI);
y_co_view = sind(VZA_ROI).*sind(VAA_ROI);
z_co_view = cosd(VZA_ROI); ... Gives height info ... Not included in BRDF

%plotting the solar angle and View Angle
figure
scatter3(x_co_solar,y_co_solar,z_co_solar)
hold on
scatter3(x_co_view,y_co_view,z_co_view)


meanSZA_ROI = nanmean(SZA_ROI);
meanSAA_ROI = nanmean(SAA_ROI);
meanVZA_ROI = nanmean(VZA_ROI);
meanVAA_ROI = nanmean(VAA_ROI);

%preallocation of the variables
% forgot to change the row vector to column vector in mean of ROI
[row,col] = size(filteredMean);
predicted = zeros(col,row);
corrected = zeros(col,row);

% performing the normality Test of the TOA Mean Reflectance
% creating the Normality Test Structure 
% checking if reflectance shows normal behaviour
%intializing the structure array
normality_Test(bands) = struct();
for band = 1:bands
    alpha  = 0.05;
    [h,p_value]  = adtest(mean_value(:,band));
    % creating structure of the normality Test of the Reflectance 
    normality_Test(band).band = strcat('Band',num2str(band));
    normality_Test(band).unfilteredDataIsNormal = logical(~h);
    normality_Test(band).unfilteredP_value = p_value;
    [h,p_value]  = adtest(filteredMean(:,band));
    normality_Test(band).filteredIsNormal= logical(~h);
    normality_Test(band).filteredP_value = p_value;
end


for band = 1:bands
    % creating the table 
    table1 = table(filteredMean(:,band),transpose(y_co_solar),...
    transpose(x_co_solar),transpose(y_co_view),...
    transpose(x_co_view));
    table1.Properties.VariableNames = {'TOA_REFLECTANCE','sy','sx','vy','vx'};
     %{
        MODEL DEFINITION
        TOA_REFLECTANCE = \beta1 + \beta2*sy + \beta3*sx + \beta4*vy +
                                 \beta5*vx
     %}
    model = fitlm(table1,'TOA_REFLECTANCE~ sy + sx + vy + vx ');
    disp(model)
    %{
    plotResiduals(model,plotype,Name,Value)
    plottype — Plot type
    'histogram' (default) | 'caseorder' | 'fitted' | 'lagged' | 'probability' | 'symmetry'
    %}
    figure 
    plotResiduals(model);
    
   %{ 
        The Residuals of all the fitted line should be normal 
    %}
    % model.Residuals is in table and kstest takes list as input
    
    % shapirowilk-test 
    alpha  = 0.05;
    % selecting the Raw data from  the residuals
    % Raw  data represents the observed minus the fitted values
    raw = 1;
    residuals = table2array(model.Residuals);
    [h,p_value,W]  = swtest(residuals(:,raw),alpha);
    if p_value > alpha
        printString =  strcat('SWTEST: The band', {' '},num2str(band),{' '},' residuals is Normal,p-value'...
            ,{' '},num2str(p_value));
        disp(printString);
    else
         printString =  strcat('SWTEST: The band',{' '}, num2str(band),{' '},' residuals is not Normal,p-value'...
            ,{' '},num2str(p_value));
        warning(printString{1});
    end
    
   % Using the Anderson Darlington  to verify normality of the Data
    alpha  = 0.05;
    [h,p_value]  = adtest(residuals(:,2));
    if ~h
        printString =  strcat('ADTEST : The band', {' '},num2str(band),{' '},' residuals is Normal,p-value'...
            ,{' '},num2str(p_value));
        disp(printString);
    else
         printString =  strcat('ADTEST :The band',{' '}, num2str(band),{' '},' residuals is not Normal,p-value'...
            ,{' '},num2str(p_value));
        warning(printString{1});
    end
    
    % Value Inflation factor 
    %vif.(band) = diag(inv(corrcoef(band_6_data)));
    
    %slope estimate 
    betahat  = model.Coefficients.Estimate;
    
    % calculation of the mean reflectance or reference reflectance
    sx  = sind(meanSZA_ROI).*cosd(meanSAA_ROI);
    sy = sind(meanSZA_ROI).*sind(meanSAA_ROI); 
    vx = sind(meanVZA_ROI).*cosd(meanVAA_ROI);
    vy = sind(meanVZA_ROI).*sind(meanVAA_ROI);
    %reference refelctance
    reference = betahat(1) + betahat(2)*sy + betahat(3)*sx + betahat(4)*vy + ...
        betahat(5)*vx;
    predicted(band,:)= betahat(1) + betahat(2)*y_co_solar + betahat(3)*x_co_solar ...
                           + betahat(4)*y_co_view + betahat(5)*x_co_view;
                       
   % transposing the filteredMena
   % FilteredMean is the column matrix  and others are row matrix
    corrected(band,:) = (transpose(filteredMean(:,band))*reference)./predicted(band,:);
end 
% transposing the corrected Reflectance to the column matrixs
correctedReflectance = transpose(corrected);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              UNCERTAINITY CALCULATION AFTER BRDF CORRECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
UNCERTAINITY is an interval around the measurement in which repeated
measurement will fail
ERROR is the degree to which a Measurement agrees with the true Value.
Error is almost never  what we are interested in. In science we typically
do noot know the 'True' Value.

Rather we are interested in the Uncertainity. This is what we need to
quantify in the measurement.

The Uncertainity, in the lab, describes the distance from your measurement
result within which your setup has determined that the true value is likely
to lie.

"Full  Width at half Maximum FWHM" another way of giving the uncetainity.
FWHM  = 2.35 \sigma
rms = \sigma

%}
% calculating the mean and std for comparision only
% not used in the calculation
meanOfFilteredMean = mean(filteredMean);
% stdOfFilteredMean = std(meanOfFilteredMean);
stdOfFilteredMean = std(filteredMean);
meanOfBrdfCorrectedTOA = mean(correctedReflectance); 
stdOfBrdfCorrection  = std(correctedReflectance);

temporalUncertainity = (stdOfMean(:,1:7)./meanOfMean(:,1:7)).*100;
filteredTemporalUncertainity  = (stdOfFilteredMean(:,1:7)./meanOfFilteredMean(:,1:7))*100;
brdfTemporalUncertainity = (stdOfBrdfCorrection(:,1:7)./meanOfBrdfCorrectedTOA(:,1:7))*100;

uncertainityName = {'Temporal Uncertainity';strcat('Filtered Temporal Uncertainity');...
                'BRDF Corrected Uncertainity'};
bands = [temporalUncertainity;filteredTemporalUncertainity;...
            brdfTemporalUncertainity];
uncertainityTable = table(uncertainityName,bands);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Title text for all the filtered plots
titleText = strcat('TOA Reflectance of Libya 4 (Landsat 8) after filtering');

% outlier removed temporal data plotting
figure
for band = 1:bands
   e = plot(filteredMean(:,9),filteredMean(:,band));
   e.Marker = 'o';
   e.MarkerSize = 10;
   e.Color = bandColors{band};
   e.LineWidth = 0.75;
   e.MarkerFaceColor = bandColors{band};
   e.MarkerEdgeColor = bandColors{band};
   % removes the connecting line between the points
   e.LineStyle = 'None';
   xlabel('Year')
   ylabel('TOA Reflectance','FontSize',16,'FontWeight','bold')
   title(titleText,'FontSize',16,'FontWeight','bold')
   legend(bandName{1:7},'location','northeast')
   hold on
   grid on
   % limit the Y-axis 
   %ylim ([0, 1])
end

%plotting the temporal data with spatial uncertainiy.
figure 
for band = 1:bands
   e = errorbar(filteredMean(:,9),filteredMean(:,band),filteredStd(:,band));
   e.Marker = 'o';
   e.MarkerSize = 10;
   e.Color = bandColors{band};
   e.LineWidth = 0.75;
   e.CapSize = 10;
   e.LineStyle = 'None';
   e.MarkerFaceColor = bandColors{band};
   e.MarkerEdgeColor = bandColors{band};
   xlabel('Year')
   ylabel('TOA Reflectance','FontSize',16,'FontWeight','bold')
   title(titleText,'FontSize',16,'FontWeight','bold')
   legend(bandName{1:7},'location','northeast')
   hold on
   grid on
end


% BRDF corrected reflectance
titleText = strcat('TOA Reflectance of Libya 4 (Landsat 8) after filtering and BRDF Correction');
figure
for band = 1:7
   e = plot(filteredMean(:,9),correctedReflectance(:,band));
   e.Marker = 'o';
   e.MarkerSize = 10;
   e.Color = bandColors{band};
   e.LineWidth = 0.75;
   e.MarkerFaceColor = bandColors{band};
   e.MarkerEdgeColor = bandColors{band};
   % removes the connecting line between the points
   e.LineStyle = 'None';
   xlabel('Year')
   ylabel('TOA Reflectance','FontSize',16,'FontWeight','bold')
   title(titleText,'FontSize',16,'FontWeight','bold')
   legend(bandName{1:7},'location','northeast')
   hold on
   grid on
   % limit the Y-axis 
   %ylim ([0, 1])
end

%plotting the temporal data with spatial uncertainiy.
figure 
for band = 1:7
   e = errorbar(filteredMean(:,9),correctedReflectance(:,band),filteredStd(:,band));
   e.Marker = 'o';
   e.MarkerSize = 10;
   e.Color = bandColors{band};
   e.LineWidth = 0.75;
   e.CapSize = 10;
   e.LineStyle = 'None';
   e.MarkerFaceColor = bandColors{band};
   e.MarkerEdgeColor = bandColors{band};
   xlabel('Year')
   ylabel('TOA Reflectance','FontSize',16,'FontWeight','bold')
   title(titleText,'FontSize',16,'FontWeight','bold')
   legend(bandName{1:7},'location','northeast')
   hold on
   grid on
end







