% function [Check] = FoodWaste()
clc;
clear all;

%% Steve revising in Nov 2017
%% Steve further revised sensitivity analysis (Fig. 2) in Feb 2018
%% Steve further revised to differentiate developed/developing countries in May 2018
%% Erin revised to add new price elasticities May 2018
%% Margaret began revising in Jan 2021 - created supply-demand curves and fixed code errors 
%% Working Code - December 2021
%% Working Code - February 2022

n_Elast_Supply = zeros(1000,7,8);
n_Elast_Demand = zeros(1000,7,8);
ConsumptionIncrease = zeros(121,7,8);
wasteAvoided = zeros(121,7,8);
lossAvoided = zeros(121,7,8);
ReboundPct = zeros(121,7,8);
Elasticities_S = zeros(121,7,8);
Elasticities_D = zeros(121,7,8);

PrintFigs = 1;
        

%% Initialize ranges of elasticities and P/Q at present for each ROI (region of interest)

% Waste data by category, phase and region (also loads arrays of names)
load('Data [Food Waste].mat')

%carbonImpact (8x4 double)

%demandElast (8x3x7 double)    
    %1D: SDG region (see regionSDG description)  
    %2D: Low (column 1), Average (column 2), High (column 3) demand
    %Elasticity estimates from Green et al (2013) applied to relevant food
    %types. See SI for additional information on how we equated Green et
    %al's food types with food types from the FAO 
    %3D: Food Type (see foodType description)
    
%foodKey (7x1 cell)
    %1D: Food Type (see foodType description)
    %2D: Food type corresponding to Green et al (2013) food types. See SI
    %for additional information. 

%foodKey2 (7x1 cell)

%foodLoss (8x1 vector)
    %1D: SDG Region (see regionSDG description)
    %2D: Food Loss value from Post-Harvest to Distribution along the food
    %supply chain in percent (%) from the 2019 FAO State of Food and
    %Agriculture (SOFA) Report Figure 3
    
%foodType (7x1 cell array)
    %Food Types: 1. Cereals, 2. Fruits and Vegetables, 3. Meat, 4,
    %Milk, 5. Oilcrops and Pulses, 6. Roots and Tubers, 7. Fish and Seafood
    
%foodWaste (8x1 cell array)
    %1D: SDG Region (see regionSDG description)
    %2D: Food waste in percent (%) of supply wasted for food
    %service, retail, and household waste from the 2021 Food Waste Index
    %(FWI) report Level 1 Annex. 
    
%landImpact (8x4 double)

%n
    
%price (8x1 vector) 
    %1D: SDG Region (See regionSDG description)
    %2D: FAO Food Price Indices across different regions from FAOSTAT
 
%regionalElasticities (8x3x3 double)
    %1D: Food Type from Green et al (2013). See foodTypeKey for how we
    %convert between Green et al food types and FAO food types
    %2D: Food demand elasticities calculated by green et al. Columns 1 and
    %3 refer to the high and low end confidence intervals of the average
    %value providd in column 2 for a given food typ
    %3D: Income level. 1. Low Income, 2. Middle Income, 3. High Income.
    %See regionSDG and regionIncomeKey for more details. 

%regionIncomeKey
    %1D: SDG region
    %2D: Income level of SDG region in high (3), middle (2), or low (1)
    %income. 
    
%regionSDG (8x2) vector
    %1D: Region Names (1. Australia and New Zealand, 2. Central and Southern
    %Asia, 3. Eastern and South-Eastern Asia, 4. Latin America and the
    %Caribbean,5. Northern America and Europe, 6. Oceania (excluding
    %Australia and New Zealand), 7. Sub-Saharan Africa, 8. Western Asia and
    %Northern Africa)
    %2D: Income classification based on World Bank and Green et al (2013).
    %Regions can be classified as either High, Middle, or Low income.
    
%supply (8x7 double) 
    %1D: SDG Region (see regionSDG description)
    %2D: Food Type (see foodType description)
    %Regional supply quantities in Megatonnes per year (Mt/yr) for each
    %region and food type combination using FAOSTAT data on production,
    %imports, and exports. 
    %Supply = Production - Exports + Imports - Losses
    
%supplyElast (8x7 double)
    %1D: SDG Region (see regionSDG description)
    %2D: Food Type (see foodType description)
    %Supply elasticities for each region and food type combination using
    %data from USDA SWOPSIM model estimates from J.V. Stout (1991) Appendix 
    %Table 7. See SI for how we equated Stout's food types with FAO food
    %types. 

%waterImpact (8x4 double)


%% Calculate region- and food-specific supply and demand curves and rebound due to waste avoided
    
for FTI=1:7 %Food Type Index (FTI) test with :1 just cereals, change to :7 to include all food types
    
    Max_Es_in = max(supplyElast(1:8,FTI)); %sets max Es, 0.35 for cereals
    Min_Es_in = min(supplyElast(1:8,FTI)); %sets min Es, 0.29 for cereals

    for SDGI=1:8 %SDG Index (SDGI) Region
        
        incomeLevel = regionIncomeKey(SDGI); %creates a 
        foodLevel = foodKey(FTI);
        
        Max_Ed_in = regionalElasticities(foodLevel,3,incomeLevel); 
        Mean_Ed_in = regionalElasticities(foodLevel,2,incomeLevel);
        Min_Ed_in = regionalElasticities(foodLevel,1,incomeLevel);
        
        %for cereals, takes min_Ed of -0.32 from Swopsim...
        Q_in(FTI,SDGI) = supply(SDGI,FTI); %Creates a vector of supply 
        P_in(FTI,SDGI) = price(FTI,1); %Global commodity price values
        
        % Repeatedly (n times) calculate additional consumption per unit of avoided waste (i.e. ReboundRatio) for the specific foodtype and ROI
        rng('default'); %for reproducibility
        counter=1;
        for i = 1:n
            %randomly draw from trinagular distribution of food- and region-specific demand elasticities
            pd = makedist('Triangular','a',Min_Ed_in,'b',Mean_Ed_in,'c',Max_Ed_in); %probability distribution
            n_Elast_Demand(counter,FTI,SDGI) = random(pd,1,1);
        
            %LHSU randomly draws from food- and region-specific range assuming uniform distribution
            n_Elast_Supply(counter,FTI,SDGI) = lhsu(Min_Es_in,Max_Es_in,1);
           
            counter=counter+1;
       
        end %for i
              
    end %for ROI--Closes the 'for' loop
    
end %for foodtype


        

%% Figure 1 (using the n-loop results...)
counter1 = 1; 

 for FTI=1:7 %Food Type Index (FDI); change to :7 to include all food types

     for SDGI= 1:8 %SDG Index (SDGI); change to :8 to include all regions
         figure % Plot shift in supply curve (i.e. Wav) and new equilibrium (including percentile range)
         hold on
         
         % Indicate desired percentiles
         percentiles = (0:10:100); %length = 11
         
         for loop1 = 1:length(percentiles)
             pct1 = percentiles(loop1);
             pct_elast_d = -abs(prctile(n_Elast_Demand(:,FTI,SDGI),pct1));

             for loop = 1:length(percentiles)
                 pct = percentiles(loop);
                 
                 %Using percentiles of the permutations generated by the "n-loop" above
                 pct_wastepct = 0.50;
                 pct_elast_s = abs(prctile(n_Elast_Supply(:,FTI,SDGI),pct));

                 %calculate supply curve
                 Ps = (P_in(FTI,SDGI)); 
                 Qs = (Q_in(FTI,SDGI));
                 Cs = Qs(1)/ (Ps(1) .^ pct_elast_s); %constant in power law function for supply curve

                 %Ps and Qs reflect initial equilibrium - from FAO index data
                 Ps(2) = (Ps(1) .* 1.1); %Make another point based on current equilibrium; change by increase 10%
                 Qs(2) = Cs * (Ps(2) ^ pct_elast_s); %Calculates Qs(2) based on a 10% increase in price using the constant C and the elasticity
              
                 Pd = (P_in(FTI,SDGI)); 
                 Qd = (Q_in(FTI,SDGI));

                 Cd = Qd(1) / (Pd(1).^(pct_elast_d)); %constant in power law function for demand curve

                 Pd(2) = (Pd(1) .* 1.1); %made another point based on current equilibrium; change in 10%
                 Qd(2) =  Cd * (Pd(2) ^ pct_elast_d);

                 %initialize range; Suppy and Demand are equal by assumption at the equilibrium
                 if Qs < 5
                    minQ = 0.01;
                    diff = 0.01;
                 else
                    minQ = 1;
                    diff = 1;
                 end 
                 
                 maxQ = max(max(Qs),max(Qd))*3;
                 D=[];
                 S=[];

                 %calculate supply and demand curves
                 [Eq_Q,Eq_P] = polyxpoly(Qs,Ps,Qd,Pd); %this checks that Eq_Q and Eq_P are the initial eq. values (from data) and renames for ease of use below

                 D(1,:) = (minQ:diff:maxQ); %initializes quantities (X values) for demand curve
                 S(1,:) = (minQ:diff:maxQ);  %initializes quantities (X values) for the supply curve
                 D(2,:) = (D(1,:)./Cd).^(1/pct_elast_d); %calculates prices based on quantities D(1,:)
                 S(2,:) = (S(1,:)./Cs).^(1/pct_elast_s); %calculates prices based on quantities S(1,:)

                 %plot present-day supply and demand curves, including calculated equilibrium

                 plot(S(1,:),S(2,:)) %Original Supply Curve
                 hold on
                 ylim([0 200])
                 scatter(Eq_Q,Eq_P,'og') 
                 text(Eq_Q,Eq_P,strcat('(',num2str(Eq_Q),',',num2str(Eq_P),')'))

                 %calculate new equilibrium quantity based on waste avoided and rebound
                 S_adj = foodLoss(SDGI,1)/100 * production(SDGI,FTI) * pct_wastepct; %adjusted by production
                 S2 = S(1,:) + S_adj; %Adjusted supply curve
                 plot(S2,S(2,:)) %plots adjusted supply curve
                
                 D_adj = foodWaste(SDGI,1)/100*Q_in(FTI,SDGI) * pct_wastepct; %do total supply
                 D2 = D(1,:) - D_adj; %Adjusted demand curve
                 plot(D2,D(2,:)) %plots adjusted demand curve
                
                 %plot adjusted supply curve, including new calculated equilibrium
                 [Eq_Q_new,Eq_P_new] = polyxpoly(S2,S(2,:),D2,D(2,:)); %New Equilibrium Price and Quantity for adjusted supply curve
                 %scatter(Eq_Q_new,Eq_P_new,'og') 
                 %text(Eq_Q_new,Eq_P_new,strcat('(',num2str(Eq_Q_new),',',num2str(Eq_P_new),')'))
                 
                 %These are the results of each foodtype, ROI and "n-loop" percentile combination

                 ConsumptionIncrease(counter1,FTI,SDGI) = Eq_Q_new - Eq_Q;
                 wasteAvoided(counter1,FTI,SDGI) = D_adj; %equivalent to the adjustment, Eq_Q_adj
                 lossAvoided(counter1,FTI,SDGI) = S_adj; 
                 ReboundPct(counter1,FTI,SDGI) = (wasteAvoided(counter1, FTI,SDGI) + ConsumptionIncrease(counter1,FTI,SDGI)) / (wasteAvoided(counter1, FTI,SDGI) + lossAvoided(counter1,FTI,SDGI));
                 Elasticities_S(counter1,FTI,SDGI) = pct_elast_s;
                 Elasticities_D(counter1,FTI,SDGI) = pct_elast_d;

                 plotname = strcat(regionSDG(SDGI),foodType(FTI)); %names the graphs with supply and demand curves
                 title(plotname)
                 xlabel('Quantity');
                 ylabel('Price');
                 counter1 = counter1 + 1; 
             
             end %percentile for loop for supply
             plot(D(1,:),D(2,:)) %Demand Curve
             
             
         end %percentile for loop for demand
         hold off
         
         %Names each graph and saves them 
         counter1 = 1; 
         savename = strcat(regionSDG(SDGI),foodType(FTI));
       
         if PrintFigs == 1
           print(gcf,'-depsc','-painters', savename{1});
         end 

    end %for ROI--Closes the 'for' loop
    
end %for foodtype


%Summary Statistics

for f = 1:7
    for r = 1:8
        Y_r(r,f) = median(ReboundPct(:,f,r));
        Y_w(r,f) = median(wasteAvoided(:,f,r));
        Y_l(r,f) = median(lossAvoided(:,f,r));
        Y_c(r,f) = median(ConsumptionIncrease(:,f,r));
        
        P_r_l(r,f) = prctile(ReboundPct(:,f,r),2.5); %Lower percentile
        P_r_u(r,f) = prctile(ReboundPct(:,f,r),97.5); %Upper percentile
        
        P_c_l(r,f) = prctile(ConsumptionIncrease(:,f,r),2.5);
        P_c_u(r,f) = prctile(ConsumptionIncrease(:,f,r),97.5);
    end
end

%Back-of-the-Envelope Emissions Calculation

for f = 1:7 
    for r = 1:8
        
        foodLevel2 = foodKey2(f);
        possibleEmissions(r,f) = (Y_w(r,f) +Y_l(r,f)) * carbonImpact(r,foodLevel2); %(Mt CO2 eqv) %Impact factor =  tonne CO2 eqv / tonne of food lost possible emissions avoided
        possibleWater(r,f) = (Y_w(r,f) +Y_l(r,f)) * 10^6 * waterImpact(r,foodLevel2); %(m^3 water) Impact Factor = m^3 water/ tonne of food
        possibleLand(r,f) = (Y_w(r,f) +Y_l(r,f)) * 10^6 * landImpact(r,foodLevel2); % (ha of land) Impact factor = ha land/ tonne of food
        
        actualEmissions(r,f) = possibleEmissions(r,f) - (Y_w(r,f) + Y_c(r,f)) * carbonImpact(r,foodLevel2); %actual emissions avoided with rebound effect
        actualWater(r,f) = possibleWater(r,f) - (Y_w(r,f) + Y_c(r,f)) * waterImpact(r,foodLevel2)*10^6;
        actualLand(r,f) = possibleLand(r,f) - (Y_w(r,f) + Y_c(r,f))* landImpact(r,foodLevel2)* 10^6;
    end
end

totalPossibleWaste = sum(Y_w,"all") + sum(Y_l,"all");
totalActualWaste = totalPossibleWaste -sum(Y_w,"all") - sum(Y_c,"all");
offsetWaste = (totalActualWaste - totalPossibleWaste) / totalPossibleWaste;

totalPossibleEmissions = sum(possibleEmissions, "all")/10^3; %(Gt CO2 eqv) Total possible emissions that could be avoided from reducing FWL
totalActualEmissions =  sum(actualEmissions, "all")/10^3; %(Gt CO2 eqv) Total actual emissions reduced with rebound effects
offsetEmissions = (totalActualEmissions - totalPossibleEmissions) / totalPossibleEmissions;

totalPossibleWater = sum(possibleWater, "all");
totalActualWater =  sum(actualWater, "all");
offsetWater = (totalActualWater - totalPossibleWater) / totalPossibleWater;

totalPossibleLand = sum(possibleLand, "all");
totalActualLand =  sum(actualLand,"all");
offsetLand = (totalActualLand - totalPossibleLand) / totalPossibleLand;

globalEmissions = 52.3; %(Gt CO2 eqv) total global emissions
foodSystemEmissions = 0.26 * globalEmissions; %(Gt CO2 eqv) food system emissions
foodWasteEmissions = 0.06 * globalEmissions; %(Gt CO2 eqv) food waste emissions
possibleReduction = totalPossibleEmissions/ foodWasteEmissions * 100; % (%) ~50%
actualReduction = totalActualEmissions/foodWasteEmissions *100;% (%)

savingsAg = totalPossibleEmissions / foodSystemEmissions * 100; 
