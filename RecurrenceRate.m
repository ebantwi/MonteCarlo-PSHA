function [ZoneBeta ZoneNo ZoneLambda ] = RecurrenceRate(Mw,Years,MinMw,MaxMw)
        [NumEntry,NumZones] = size(Mw);
        PeriodforZone = length(Years);
        rate = zeros(size(Mw));
        ZoneBeta = zeros(1,NumZones)
        ZoneNo = zeros(1,NumZones)
        ZoneLambda = zeros(1,NumZones)
%          CellTruncatedMW = struct;
        for i = 1:NumZones
            CurrentMw = Mw(:,i);
            TruncatedMW  = CurrentMw((Mw(:,i) > 0))
            for j = 1:length(TruncatedMW)
                rate(j,i) = sum(TruncatedMW >= TruncatedMW(j))./Years(i)
            end
            RateZone = rate(1:length(TruncatedMW),i)
            CurrentRate = RateZone;
            logRateZone = log10(RateZone);
            figure(i)
            semilogy(TruncatedMW,logRateZone,'ro')
            grid on
            a = fitlm(TruncatedMW , logRateZone);
            Intercept = table2array(a.Coefficients(1,1));
            Slope =  table2array(a.Coefficients(2,1));
            No = 10^Intercept
            Lambda = 10^(Intercept+Slope*MinMw)
            t = Lambda/No
            
            beta = 0.0001;
           
            while t -((exp(-beta*MinMw)-exp(-beta*MaxMw))/(1-exp(-beta*MaxMw))) < 0
                beta = beta + 0.00001
            end
            ZoneBeta(1,i) = beta;
            ZoneNo(1,i) = No;
            ZoneLambda(1,i) = Lambda;
            
%             if isfield(CellTruncatedMW, 'CurrentRate')
%             CellTruncatedMW.CurrentRate = vertcat(CellTruncatedMW.CurrentRate,mat2cell(CurrentRate,size(CurrentRate)))
%             else
%                 CellTruncatedMW.CurrentRate = mat2cell(RateZone)
%             end
        end
        
    end