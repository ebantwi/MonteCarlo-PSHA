function [ LogPSA ] = GMPEEstimates(Mag,Location,Site,GMPEConsidered,GMPEWeights)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
ValidLocationIndex  = find(Location(:,1))

%finding coordinates of epicenters within zone boundaries
LocationX = Location(:,1);
ValidLocationX = LocationX(ValidLocationIndex)
LocationY = Location(:,2);
ValidLocationY = LocationY(ValidLocationIndex)

%coordinates of epicenters within boundaries
ValidLocation = horzcat(ValidLocationX,ValidLocationY);

%corresponding magnitudes of valid epicenters
ValidMag = Mag(ValidLocationIndex')
wgs84 = wgs84Ellipsoid('meters');
SitelatBoundary = Site(1,1);
SitelonBoundary = Site(1,2);
[xSiteLocation,ySiteLocation] = geodetic2ecef(wgs84,SitelatBoundary,SitelonBoundary,200);

% coordinates of site of interest
Site = [xSiteLocation,ySiteLocation];

%initializing source to site distance as a column vector
SourceToSiteDistance = zeros(length(ValidLocation),1)

%distance from epicenters to site of interest
for i = 1:length(ValidLocation)
    SourceToSiteDistance(i,1) = (sum((ValidLocation(i,:) - Site).^2).^0.5)/1000 %in kilometres
end

[SizeGMPEConsidered,~] = size(GMPEConsidered); %This should be in a column vector
[SizeGMPEWeights,~] = size(GMPEWeights); %This should be in a column vector

load GMPE_data.mat
for j = 1:SizeGMPEConsidered 
    Name = cell2mat(GMPEConsidered(j,:))
    Weight = GMPEWeights(j,1);
   
    switch Name

        case 'AtkinsonHardrock'
           [PeriodInterested Coefficients] = size(AtkinsonHardrock)
           AtkinsonHardrock = horzcat(AtkinsonHardrock,ones(PeriodInterested,1));
           AtkinsonHardrockFinal = AtkinsonHardrock(:,2:end);
           MeanGMPEHR = zeros(length(ValidLocation),1)
           
           for l = 1:length(ValidLocation)
               % starting from one, select for loop all magnitudes of valid
               % epicenters
               M = ValidMag(l,1)
               %starting from one, select for loop all distances tp various
               %epicenters from site.
               Rcd = SourceToSiteDistance(l,1)
               
               f0 = max(log(10/Rcd),0);
               f1 = min(log(Rcd),log(70));
               f2 = max(log(Rcd),140);
               S = 0
              Vector = [1 M M^2 f1 M*f1 f2 M*f2 f0 M*f0 Rcd S];
              if MeanGMPEHR == 0
              MeanGMPEHR = Vector*(AtkinsonHardrockFinal')*Weight; 
              else
                   CurrentMeanGMPEHR = Vector*(AtkinsonHardrockFinal')*Weight; 
                   MeanGMPEHR = vertcat(MeanGMPEHR,CurrentMeanGMPEHR);
              end    
           
           end
           
        case 'AtkinsonBC'
           [PeriodInterested Coefficients] = size(AtkinsonBC)
           AtkinsonBC = horzcat(AtkinsonBC,ones(PeriodInterested,1));
           AtkinsonBCFinal = AtkinsonBC(:,2:end);
           MeanGMPEBC = zeros(length(ValidLocation),1)
         
           for l = 1:length(ValidLocation)
               M = ValidMag(l,1)
               Rcd = SourceToSiteDistance(l,1)
               f0 = max(log(10/Rcd),0);
               f1 = min(log(Rcd),log(70));
               f2 = max(log(Rcd),140);
               S = 0
              Vector = [1 M M^2 f1 M*f1 f2 M*f2 f0 M*f0 Rcd S];
              if MeanGMPEBC == 0
              MeanGMPEBC =Vector*(AtkinsonBCFinal')*Weight; 
              else
                   CurrentMeanGMPEBC = Vector*(AtkinsonBCFinal')*Weight; 
                   MeanGMPEBC = vertcat(MeanGMPEBC,CurrentMeanGMPEBC);
              end    
           
           end
           
       case 'PezeshkHybridEM'
   %line below gives the number of periods used and number of coefficients available 
          [PeriodInterested Coefficients] = size(PezeshkHybridEM)
   %next line is to pick data in all rows from c1-c10 for Pezeshk GMPE computation        
          PezeshkHybridEMFinal = PezeshkHybridEM(:,2:11);
    % preallocate MeanGMPEPezeshk with a matrix of size length(Valid locations),1      
           MeanGMPEPezeshk = zeros(length(ValidLocation),1) 
           c11 = PezeshkHybridEM(:,12)
          % matrix property: square of a matrix = (transpose of matrix) * matrix
           c11squared = (c11')*c11;
             
         for l = 1:length(ValidLocation)
               M = ValidMag(l,1) %validMag is a column vector
               
               Rrup = SourceToSiteDistance(l,1) %Rrup is a column vector
               R = ((sqrt(Rrup)^2)+ c11squared)^0.5;
               f0 = min(log(R),log(70));                        
               f1 = max(min(log(R/70),log(140/70)),0);
               f2 = max(log(R/140),0);
               
              Vector = [1 M M^2 f0 M*f0 f1 M*f1 f2 M*f2 R];
              if MeanGMPEPezeshk == 0
             MeanGMPEPezeshk = Vector*(PezeshkHybridEMFinal')*Weight; 
              else
                   CurrentMeanGMPEPezeshk = Vector*(PezeshkHybridEMFinal')*Weight; 
                   MeanGMPEPezeshk = vertcat(MeanGMPEPezeshk,CurrentMeanGMPEPezeshk);
              end    
           
           end
                  
    case 'SilvaDCM'  
   %line below gives the number of periods used and number of coefficients available 
          [PeriodInterested Coefficients] = size(SilvaDCM)
   %next line is to pick data in all rows from c1-c10 for Silva Double Corner Model GMPE computation        
%           SilvaDCMFinal = SilvaDCM(:,2:end);
%    % preallocate MeanGMPEPezeshk with a matrix of size length(Valid locations),1      
%            MeanGMPESilvaDCM = zeros(length(ValidLocation),1) 
           % log(y) = c1 +c2M + (c6+c7M)*log(R + e^c4) + c10(M-6)
           % a = (c1+c2M)  b= (c6+c7M) c = log(R+e^c4) d = c10(M-6)
           % log(y) = a + (b.*c) + d
% load('GMPE_data.mat')
          %call variables from GMPE data table of double corner model
           c1d = SilvaDCM(:,2);
           c2d = SilvaDCM(:,3);
           c4d = SilvaDCM(:,4);
           c5d = SilvaDCM(:,5);   
           c6d = SilvaDCM(:,6);                   
           c7d = SilvaDCM(:,7);
           c8d = SilvaDCM(:,8);
           c10d = SilvaDCM(:,9);

% for i = 1:length(mattry)
%     a = mattry(i,2)*2 +mattry(i,3)*1 ;
%      b = sum(a(:,1:end));
%  c(1,i) = b

  MeanGMPESilvaDCM = zeros(length(ValidLocation),1) 

for l = 1: length(ValidLocation)
    M = ValidMag(l,1) 
    R = SourceToSiteDistance(l,1)
    a = c1d + M*c2d;
    b = c6d + M*c7d;
    c = log(R + exp(c4d));
    d = ((M-6)^2)*c10d;
    e = b.*c
    if  MeanGMPESilvaDCM ==0
     MeanGMPESilvaDCM = (a+ (b.*c) +d)'*Weight;
   else
  CurrentMeanGMPESilvaDCM = (a+ (b.*c) +d)'*Weight;
  MeanGMPESilvaDCM = vertcat(MeanGMPESilvaDCM, CurrentMeanGMPESilvaDCM);         
    end
            
end


        case 'SilvaDCMSaturation'
   %line below gives the number of periods used and number of coefficients available 
          [PeriodInterested Coefficients] = size(SilvaDCMSaturation)
   %next line is to pick data in all rows from c1-c10 for Silva Double Corner Model GMPE computation with saturation       
%           SilvaDCMFinal = SilvaDCM(:,2:end);
%    % preallocate MeanGMPEPezeshk with a matrix of size length(Valid locations),1      
%            MeanGMPESilvaDCM = zeros(length(ValidLocation),1) 
           % log(y) = c1 +c2M + (c6+c7M)*log(R + e^c4) + c10(M-6)
           % a = (c1+c2M)  b= (c6+c7M) c = log(R+e^c4) d = c10(M-6)
           % log(y) = a + (b.*c) + d
% load('GMPE_data.mat')
          %call variables from GMPE data table of double corner model
           c1d = SilvaDCMSaturation(:,2);
           c2d = SilvaDCMSaturation(:,3);
           c4d = SilvaDCMSaturation(:,4);
           c5d = SilvaDCMSaturation(:,5);   
           c6d = SilvaDCMSaturation(:,6);                   
           c7d = SilvaDCMSaturation(:,7);
           c8d = SilvaDCMSaturation(:,8);
           c10d =SilvaDCMSaturation(:,9);

% for i = 1:length(mattry)
%     a = mattry(i,2)*2 +mattry(i,3)*1 ;
%      b = sum(a(:,1:end));
%  c(1,i) = b

  MeanGMPESilvaDCMSaturation = zeros(length(ValidLocation),1) 

for l = 1: length(ValidLocation)
    M = ValidMag(l,1) 
    R = SourceToSiteDistance(l,1)
    a = c1d + M*c2d;
    b = c6d + M*c7d;
    c = log(R + exp(c4d));
    d = ((M-6)^2)*c10d;
    if  MeanGMPESilvaDCMSaturation ==0
     MeanGMPESilvaDCMSaturation = (a+ (b.*c) +d)'*Weight;
   else
  CurrentMeanGMPESilvaDCMSaturation = (a+ (b.*c) +d)'*Weight;
  MeanGMPESilvaDCMSaturation = vertcat(MeanGMPESilvaDCMSaturation, CurrentMeanGMPESilvaDCMSaturation);         
    end
end

        case 'SilvaScmConStressSatu'
   %line below gives the number of periods used and number of coefficients available 
          [PeriodInterested Coefficients] = size(SilvaScmConStressSatu)
   %next line is to pick data in all rows from c1-c10 for Silva Double Corner Model GMPE computation        
%           SilvaDCMFinal = SilvaDCM(:,2:end);
%    % preallocate MeanGMPEPezeshk with a matrix of size length(Valid locations),1      
%            MeanGMPESilvaDCM = zeros(length(ValidLocation),1) 
           % log(y) = c1 +c2M + (c6+c7M)*log(R + e^c4) + c10(M-6)
           % a = (c1+c2M)  b= (c6+c7M) c = log(R+e^c4) d = c10(M-6)
           % log(y) = a + (b.*c) + d
% load('GMPE_data.mat')
          %call variables from GMPE data table of double corner model
           c1d = SilvaScmConStressSatu(:,2);
           c2d = SilvaScmConStressSatu(:,3);
           c4d = SilvaScmConStressSatu(:,4);
           c5d = SilvaScmConStressSatu(:,5);   
           c6d = SilvaScmConStressSatu(:,6);                   
           c7d = SilvaScmConStressSatu(:,7);
           c8d = SilvaScmConStressSatu(:,8);
           c10d =SilvaScmConStressSatu(:,9);

% for i = 1:length(mattry)
%     a = mattry(i,2)*2 +mattry(i,3)*1 ;
%      b = sum(a(:,1:end));
%  c(1,i) = b

  MeanGMPESilvaScmConStressSatu = zeros(length(ValidLocation),1) 

for l = 1: length(ValidLocation)
    M = ValidMag(l,1) 
    R = SourceToSiteDistance(l,1)
    a = c1d + M*c2d;
    b = c6d + M*c7d;
    c = log(R + exp(c4d));
    d = ((M-6)^2)*c10d;
    if  MeanGMPESilvaScmConStressSatu ==0
     MeanGMPESilvaScmConStressSatu = (a+ (b.*c) +d)'*Weight;
   else
  CurrentMeanGMPESilvaScmConStressSatu = (a+ (b.*c) +d)'*Weight;
  MeanGMPESilvaScmConStressSatu = vertcat(MeanGMPESilvaScmConStressSatu, CurrentMeanGMPESilvaScmConStressSatu);         
    end
            
end
   
    case 'SilvaScmVarStress'  
   %line below gives the number of periods used and number of coefficients available 
          [PeriodInterested Coefficients] = size(SilvaScmVarStress)
   %next line is to pick data in all rows from c1-c10 for Silva Double Corner Model GMPE computation        
%           SilvaDCMFinal = SilvaDCM(:,2:end);
%    % preallocate MeanGMPEPezeshk with a matrix of size length(Valid locations),1      
%            MeanGMPESilvaDCM = zeros(length(ValidLocation),1) 
           % log(y) = c1 +c2M + (c6+c7M)*log(R + e^c4) + c10(M-6)
           % a = (c1+c2M)  b= (c6+c7M) c = log(R+e^c4) d = c10(M-6)
           % log(y) = a + (b.*c) + d
% load('GMPE_data.mat')
          %call variables from GMPE data table of double corner model
           c1d = SilvaScmVarStress(:,2);
           c2d = SilvaScmVarStress(:,3);
           c4d = SilvaScmVarStress(:,4);
           c5d = SilvaScmVarStress(:,5);   
           c6d = SilvaScmVarStress(:,6);                   
           c7d = SilvaScmVarStress(:,7);
           c8d = SilvaScmVarStress(:,8);
           c10d =SilvaScmVarStress(:,9);

% for i = 1:length(mattry)
%     a = mattry(i,2)*2 +mattry(i,3)*1 ;
%      b = sum(a(:,1:end));
%  c(1,i) = b

  MeanGMPESilvaScmVarStress = zeros(length(ValidLocation),1) 

for l = 1: length(ValidLocation)
    M = ValidMag(l,1) 
    R = SourceToSiteDistance(l,1)
    a = c1d + M*c2d;
    b = c6d + M*c7d;
    c = log(R + exp(c4d));
    d = ((M-6)^2)*c10d;
    if  MeanGMPESilvaScmVarStress ==0
     MeanGMPESilvaScmVarStress = (a+ (b.*c) +d)'*Weight;
   else
  CurrentMeanGMPESilvaScmVarStress = (a+ (b.*c) +d)'*Weight;
  MeanGMPESilvaScmVarStress = vertcat(MeanGMPESilvaScmVarStress, CurrentMeanGMPESilvaScmVarStress);         
    end
end

        case 'SilvaScmConstStress'
   %line below gives the number of periods used and number of coefficients available 
          [PeriodInterested Coefficients] = size(SilvaScmConstStress)
   %next line is to pick data in all rows from c1-c10 for Silva Double Corner Model GMPE computation        
%           SilvaDCMFinal = SilvaDCM(:,2:end);
%    % preallocate MeanGMPEPezeshk with a matrix of size length(Valid locations),1      
%            MeanGMPESilvaDCM = zeros(length(ValidLocation),1) 
           % log(y) = c1 +c2M + (c6+c7M)*log(R + e^c4) + c10(M-6)
           % a = (c1+c2M)  b= (c6+c7M) c = log(R+e^c4) d = c10(M-6)
           % log(y) = a + (b.*c) + d
% load('GMPE_data.mat')
          %call variables from GMPE data table of double corner model
           c1d = SilvaScmConstStress(:,2);
           c2d = SilvaScmConstStress(:,3);
           c4d = SilvaScmConstStress(:,4);
           c5d = SilvaScmConstStress(:,5);   
           c6d = SilvaScmConstStress(:,6);                   
           c7d = SilvaScmConstStress(:,7);
           c8d = SilvaScmConstStress(:,8);
           c10d =SilvaScmConstStress(:,9);

% for i = 1:length(mattry)
%     a = mattry(i,2)*2 +mattry(i,3)*1 ;
%      b = sum(a(:,1:end));
%  c(1,i) = b

  MeanGMPESilvaScmConstStress = zeros(length(ValidLocation),1) 

for l = 1: length(ValidLocation)
    M = ValidMag(l,1) 
    R = SourceToSiteDistance(l,1)
    a = c1d + M*c2d;
    b = c6d + M*c7d;
    c = log(R + exp(c4d));
    d = ((M-6)^2)*c10d;
    if  MeanGMPESilvaScmConstStress ==0
     MeanGMPESilvaScmConstStress = (a+ (b.*c) +d)'*Weight;
   else
  CurrentMeanGMPESilvaScmConstStress = (a+ (b.*c) +d)'*Weight;
  MeanGMPESilvaScmConstStress = vertcat(MeanGMPESilvaScmConstStress, CurrentMeanGMPESilvaScmConstStress);         
    end
 end 

    end
end

% matrix dimensions for various GMPEs not the same.

LogPSA = MeanGMPEBC + MeanGMPEHR + MeanGMPEPezeshk + MeanGMPESilvaDCM + MeanGMPESilvaDCMSaturation + MeanGMPESilvaScmConStressSatu + MeanGMPESilvaScmVarStress + MeanGMPESilvaScmConstStress ;
