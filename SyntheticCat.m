function [CombinedTime CombinedMag CombinedLocation ] = SyntheticCat( ZoneLambda,ZoneBeta,MinMw,MaxMw,ZoneBoundary,SeismicBoundary,NumZones,SimPeriod,NumRuns)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here   
% time = zeros(SimPeriod*NumRuns,NumZones);
% Mag = zeros(SimPeriod*NumRuns,NumZones);
time = zeros(SimPeriod,NumZones);
Mag = zeros(SimPeriod,NumZones);
Currenttime = zeros(SimPeriod,NumZones);
CurrentMag = zeros(SimPeriod,NumZones);
% Location = zeros(SimPeriod*NumRuns,NumZones*2);
Location = zeros(SimPeriod,NumZones*2);
CurrentLocation = zeros(SimPeriod,NumZones*2);
for i = NumZones;
    B= ZoneBeta(:,i);
    lm = ZoneLambda(:,i);
    for j = 1:NumRuns
        if j == 1
             for x = 1:SimPeriod
        u = rand;
        t = (-log(1-u)./lm);
        m = -log(exp(-B.*MinMw)-(1-u)*(exp(-B.*MinMw)-exp(-B.*MaxMw)));
        latZone = ZoneBoundary(:,2);
        lonZone = ZoneBoundary(:,1);
        latBoundary = SeismicBoundary(:,2);
        lonBoundary = SeismicBoundary(:,1);
        wgs84 = wgs84Ellipsoid('meters');
        [xBoundary,yBoundary] = geodetic2ecef(wgs84,latBoundary,lonBoundary,200);
        [xZone,yZone]= geodetic2ecef(wgs84,latZone,lonZone,200);
        % plot(xBoundary,yBoundary,'bo')
        
        %plot(xZone,yZone,'ro')
        % ZoneShape = alphaShape(xZone,yZone,0,'HoleThreshold',1)
        k = boundary(xZone,yZone);
        %plot(xZone(k),yZone(k));
        
        kRegion = boundary(xBoundary,yBoundary);
        xBoundmin = min(xBoundary);
        xBoundmax = max(xBoundary);
        xGap = xBoundmax-xBoundmin;
        yBoundmin = min(yBoundary);
        yBoundmax = max(yBoundary);
        yGap = yBoundmax-yBoundmin;
        randpoint = rand(1,2);
        
        xValue = xBoundmin + randpoint(1,1).*xGap;
        
        yValue = yBoundmin + randpoint(1,2).*yGap;
        xyMax = [xBoundmax yBoundmax];
        xyMin = [xBoundmin yBoundmin];
        pointxy =[xValue yValue];
        
        % if inpolygon(pointxy,xBoundary(kRegion),yBoundary(kRegion))
        %     if inpolygon(pointxy,xZone(k),yZone(k))
        %         Epicenter = pointxy
        %     else
        %         Epicenter = [1,1]
        %     end
        % else
        %     Epicenter = [0,0]
        % end
        
        
        while ~(inpolygon(pointxy(1,1),pointxy(1,2),xBoundary(kRegion),yBoundary(kRegion)))
            randpoint = rand(1,2);
            
            xValue = xBoundmin + randpoint(1,1).*xGap;
            
            yValue = yBoundmin + randpoint(1,2).*yGap;
            xyMax = [xBoundmax yBoundmax];
            xyMin = [xBoundmin yBoundmin];
            pointxy =[xValue yValue];
        end
        if inpolygon(pointxy(1,1),pointxy(1,2),xZone(k),yZone(k))
            Epicenter = pointxy
        else
            Epicenter = [0,0]
        end
        
   
          


        %plot(xValue,yValue,'r*')
        %IndexMinx = find(xBoundary==(xBoundmin));
        %yBound = yBoundary(IndexMinx);
        %hold on
%         plot(xBoundary,yBoundary,'bo')
%         %plot(xBoundary(kRegion),yBoundary(kRegion));
%         Inside = inpolygon(6.35e6,0.5e5,xBoundary(kRegion),yBoundary(kRegion));

        time(x,i) = t;
    Mag(x,i) = m;
    Location(x,2*i-1) = Epicenter(1,1);
    Location(x,2*i) = Epicenter(1,2); 
    CombinedTime  = time;
    CombinedLocation  = Location;
    CombinedMag  = Mag;
    
        
        
        %SyntheticCat( [0.0888 0.1165],[0.3567 1.333],3,7,ZoneBoundary,SeismicBoundary,[1 2])
             end
        else
            for x = 1:SimPeriod
        u = rand;
        t = (-log(1-u)./lm);
        m = -log(exp(-B.*MinMw)-(1-u)*(exp(-B.*MinMw)-exp(-B.*MaxMw)));
        latZone = ZoneBoundary(:,2);
        lonZone = ZoneBoundary(:,1);
        latBoundary = SeismicBoundary(:,2);
        lonBoundary = SeismicBoundary(:,1);
        wgs84 = wgs84Ellipsoid('meters');
        [xBoundary,yBoundary] = geodetic2ecef(wgs84,latBoundary,lonBoundary,200);
        [xZone,yZone]= geodetic2ecef(wgs84,latZone,lonZone,200);
        % plot(xBoundary,yBoundary,'bo')
        
        %plot(xZone,yZone,'ro')
        % ZoneShape = alphaShape(xZone,yZone,0,'HoleThreshold',1)
        k = boundary(xZone,yZone);
        %plot(xZone(k),yZone(k));
        
        kRegion = boundary(xBoundary,yBoundary);
        xBoundmin = min(xBoundary);
        xBoundmax = max(xBoundary);
        xGap = xBoundmax-xBoundmin;
        yBoundmin = min(yBoundary);
        yBoundmax = max(yBoundary);
        yGap = yBoundmax-yBoundmin;
        randpoint = rand(1,2);
        
        xValue = xBoundmin + randpoint(1,1).*xGap;
        
        yValue = yBoundmin + randpoint(1,2).*yGap;
        xyMax = [xBoundmax yBoundmax];
        xyMin = [xBoundmin yBoundmin];
        pointxy =[xValue yValue];
        
        % if inpolygon(pointxy,xBoundary(kRegion),yBoundary(kRegion))
        %     if inpolygon(pointxy,xZone(k),yZone(k))
        %         Epicenter = pointxy
        %     else
        %         Epicenter = [1,1]
        %     end
        % else
        %     Epicenter = [0,0]
        % end
        
        
        while ~(inpolygon(pointxy(1,1),pointxy(1,2),xBoundary(kRegion),yBoundary(kRegion)))
            randpoint = rand(1,2);
            
            xValue = xBoundmin + randpoint(1,1).*xGap;
            
            yValue = yBoundmin + randpoint(1,2).*yGap;
            xyMax = [xBoundmax yBoundmax];
            xyMin = [xBoundmin yBoundmin];
            pointxy =[xValue yValue];
        end
        if inpolygon(pointxy(1,1),pointxy(1,2),xZone(k),yZone(k))
            Epicenter = pointxy
        else
            Epicenter = [0,0]
        end
        
   
          


        %plot(xValue,yValue,'r*')
        %IndexMinx = find(xBoundary==(xBoundmin));
        %yBound = yBoundary(IndexMinx);
        %hold on
%         plot(xBoundary,yBoundary,'bo')
%         %plot(xBoundary(kRegion),yBoundary(kRegion));
%         Inside = inpolygon(6.35e6,0.5e5,xBoundary(kRegion),yBoundary(kRegion));

        Currenttime(x,i) = t;
    CurrentMag(x,i) = m;
    CurrentLocation(x,2*i-1) = Epicenter(1,1);
    CurrentLocation(x,2*i) = Epicenter(1,2); 
    CombinedTime  = vertcat(time,Currenttime)
    CombinedLocation  = vertcat(Location,CurrentLocation)
    CombinedMag  = vertcat(Mag,CurrentMag)
    
        end
        
    
    
     
    end
    Mag = CombinedMag;
    Location = CombinedLocation;
    time = CombinedTime; 
   
end

end

