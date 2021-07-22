function FiducOutput = findfiducial4(x1,y1,beatn,plotf)

%%%%%%%%%%%%%%%%%%%%%

%%% smoothing LP filter

% dt = mean(diff(x1))/4;
% 
% x12=0:dt:x1(end);
% 
% y12 = interp1(x1,y1,x12);

% y12f = filtriraj(y12, 1/dt, 40, 'l');  % lowpass filter
% y12f = y12f-mean(y12f);

% hold on
% plot(x12,y12,'rx-')

x12 = x1;
y12f = y1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% max and min of the filtered signals

if max(y12f)<0
    FiducOutput.CL = 0;
    FiducOutput.dVdtmax = NaN;
    FiducOutput.dVdtmin = NaN;
    FiducOutput.PP = NaN;
    FiducOutput.MDP1 = NaN;
    FiducOutput.MDP2 = NaN;
    FiducOutput.APA = NaN;
    FiducOutput.APD90 = NaN;
    FiducOutput.APD90corr = NaN;
    FiducOutput.APD50 = NaN;
    FiducOutput.APD20 = NaN;
    FiducOutput.DDR = NaN;
    FiducOutput.DDR100 = NaN;
    
    return
end

[peaksy, peaksx_ind] = findpeaks(y12f,'minpeakheight',0);

if length(peaksx_ind)<2
        
    FiducOutput.CL = 0;
    FiducOutput.dVdtmax = NaN;
    FiducOutput.dVdtmin = NaN;
    FiducOutput.PP = NaN;
    FiducOutput.MDP1 = NaN;
    FiducOutput.MDP2 = NaN;
    FiducOutput.APA = NaN;
    FiducOutput.APD90 = NaN;
    FiducOutput.APD90corr = NaN;
    FiducOutput.APD50 = NaN;
    FiducOutput.APD20 = NaN;
    FiducOutput.DDR = NaN;
    FiducOutput.DDR100 = NaN;

    return
end

% choose the peak

if beatn == 0
    chosenbeat = length(peaksx_ind)-1;
else chosenbeat = beatn;
end

peakx = peaksx_ind(chosenbeat);
peaky = y12f(peakx); % also max AP -- to calculate APD50 etc.


%%% cycle length between the last and before the last peak

CL = x12(peaksx_ind(chosenbeat+1))-x12(peaksx_ind(chosenbeat));

%%% minima between the last three peaks

[peakmin1, peakminind1] = min(y12f(peaksx_ind(chosenbeat):peaksx_ind(chosenbeat+1)));
[peakmin_1, peakminind_1] = min(y12f(peaksx_ind(chosenbeat-1):peaksx_ind(chosenbeat)));

peakminind1 = peakminind1 + peaksx_ind(chosenbeat) -1;
peakminind_1 = peakminind_1 + peaksx_ind(chosenbeat-1) -1;

%%% find APD at different levels between the one before last peak and two
%%% minima

APA = y12f(peakx)-y12f(peakminind1);

maxAP90 = y12f(peakx)-0.9*APA;
maxAP50 = y12f(peakx)-0.5*APA;
maxAP20 = y12f(peakx)-0.2*APA;

AP90ind_back_r = find(y12f(peakminind_1:peakx)<maxAP90,1,'last');
AP90ind_back = AP90ind_back_r + peakminind_1 - 1;

AP50ind_back_r = find(y12f(peakminind_1:peakx)<maxAP50,1,'last');
AP50ind_back = AP50ind_back_r + peakminind_1 - 1;

AP20ind_back_r = find(y12f(peakminind_1:peakx)<maxAP20,1,'last');
AP20ind_back = AP20ind_back_r + peakminind_1 - 1;

%%%

AP90ind_front_r = find(y12f(peakx:peakminind1)<maxAP90,1,'first');
AP90ind_front = AP90ind_front_r + peakx - 1;

AP50ind_front_r = find(y12f(peakx:peakminind1)<maxAP50,1,'first');
AP50ind_front = AP50ind_front_r + peakx - 1;

AP20ind_front_r = find(y12f(peakx:peakminind1)<maxAP20,1,'first');
AP20ind_front = AP20ind_front_r + peakx - 1;

APD90 = x12(AP90ind_front)-x12(AP90ind_back);
APD50 = x12(AP50ind_front)-x12(AP50ind_back);
APD20 = x12(AP20ind_front)-x12(AP20ind_back);
APD90corr = x12(AP90ind_front)-x12(AP50ind_back);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotf
    f1=figure;plot(x12,y12f)
    hold on
    plot(x12(peakminind_1),y12f(peakminind_1),'rx','linewidth',2)
    plot(x12(peakminind1),y12f(peakminind1),'rx','linewidth',2)
    plot(x12(AP90ind_back),y12f(AP90ind_back),'gx','linewidth',2)
    plot(x12(AP90ind_front),y12f(AP90ind_front),'gx','linewidth',2)
    plot(x12(AP50ind_back),y12f(AP50ind_back),'gx','linewidth',2)
    plot(x12(AP50ind_front),y12f(AP50ind_front),'gx','linewidth',2)
    plot(x12(AP20ind_back),y12f(AP20ind_back),'gx','linewidth',2)
    plot(x12(AP20ind_front),y12f(AP20ind_front),'gx','linewidth',2)
    
    plot(x12(1:end-1),diff(y12f)*50)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mdp_ap50_int = x12(AP50ind_back)-x12(peakminind_1);

segment_rel = [0.4 0.8]*mdp_ap50_int;
segment_DDx = x12(peakminind_1)+segment_rel;

segment_DD100x = [x12(peakminind_1)+0.02 x12(peakminind_1)+0.1];

DDx_ind1 = find(x12<segment_DDx(1),1,'last');
DDx_ind2 = find(x12>segment_DDx(2),1,'first');

DD100x_ind1 = find(x12<segment_DD100x(1),1,'last');
DD100x_ind2 = find(x12>segment_DD100x(2),1,'first');

%%% linear regression coefficients

xvect1 = x12(DDx_ind1:DDx_ind2);
yvect1 = y12f(DDx_ind1:DDx_ind2);
xvect1c = [ones(length(xvect1),1) xvect1];

k1 = xvect1c\yvect1;
ycalc1 = xvect1c*k1;

xvect2 = x12(DD100x_ind1:DD100x_ind2);
yvect2 = y12f(DD100x_ind1:DD100x_ind2);
xvect2c = [ones(length(xvect2),1) xvect2];

k2 = xvect2c\yvect2;
ycalc2 = xvect2c*k2;

if plotf
    plot(x12(DDx_ind1:DDx_ind2),y12f(DDx_ind1:DDx_ind2),'r-','linewidth',2)
    plot(x12(DD100x_ind1:DD100x_ind2),y12f(DD100x_ind1:DD100x_ind2),'r-','linewidth',2)

    f2=figure;
    plot(xvect1,yvect1,'rx');hold on
    plot(xvect1,ycalc1)

    plot(xvect2,yvect2,'bx')
    plot(xvect2,ycalc2)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%5

maxdVdt = 0;
maxdVdt_x = [];

mindVdt = 0;
mindVdt_x = [];

for count1 = peakminind_1:peakx
    tempdV = y12f(count1+1)-y12f(count1);
    tempdt = x12(count1+1)-x12(count1);
    tempdVdt = tempdV/tempdt;
    
    if tempdVdt>maxdVdt
        maxdVdt = tempdVdt;
        maxdVdt_x = x12(count1);
        maxdVdt_count = count1;
    end
        
end

clear tempd*

for count2 = peakx:peakminind1
    tempdV = y12f(count2+1)-y12f(count2);
    tempdt = x12(count2+1)-x12(count2);
    tempdVdt = tempdV/tempdt;
        
    if tempdVdt<mindVdt
        mindVdt = tempdVdt;
        mindVdt_x = x12(count2);
        mindVdt_count = count2;
    end
    
end

if plotf
    figure(f1)
    plot(mindVdt_x,y12f(mindVdt_count),'rx')
    plot(maxdVdt_x,y12f(maxdVdt_count),'rx')
end

FiducOutput.CL = CL;
FiducOutput.dVdtmax = maxdVdt;
FiducOutput.dVdtmin = mindVdt;
FiducOutput.PP = peaky;
FiducOutput.MDP1 = peakmin_1; % minimum after the peak
FiducOutput.MDP2 = peakmin1; % minimum before the peak
FiducOutput.APA = APA; % PP-MDP2
% FiducOutput.dMDP = FiducOutput.MDP2-(-57.1055);
FiducOutput.APD90 = APD90;
FiducOutput.APD90corr = APD90corr;
FiducOutput.APD50 = APD50;
FiducOutput.APD20 = APD20;
FiducOutput.DDR =round(k1(2)*1000)/1000;
FiducOutput.DDR100 =round(k2(2)*1000)/1000;

format short

end
% FiducOutput.APD = mindVdt_x - maxdVdt_x;
    







% 
% 
% peaksx = x12(peaksx_ind);
% 
% % figure
% % plot(x12,y12f)
% % hold on
% % plot(x12(peaksx_ind),y12f(peaksx_ind),'rx')
% 
% % minima
% 
% peakdown_ind = [];
% 
% for st = 1:length(peaksx_ind)-1
%     [miny, minx] = min(y12f(peaksx_ind(st):peaksx_ind(st+1)));
%     peakdown_ind = [peakdown_ind peaksx_ind(st)+minx-1];        
% end
% clear st
% 
% % plot(x12(peakdown_ind),y12f(peakdown_ind),'rx')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% dy12f = diff(y12f);
% 
% mindy_x = [];
% mindy_y = [];
% maxdy_x = [];
% maxdy_y = [];
% 
% for st1 = 1:length(peaksx_ind)-1
%     [miny1, minx1] = min(dy12f(peaksx_ind(st1):peaksx_ind(st1+1)));
%     mindy_x = [mindy_x minx1+peaksx_ind(st1)-1];
%     mindy_y = [mindy_y miny1];
%     
%     [maxy1, maxx1] = max(dy12f(peaksx_ind(st1):peaksx_ind(st1+1)));
%     maxdy_x = [maxdy_x maxx1+peaksx_ind(st1)-1];
%     maxdy_y = [maxdy_y maxy1];
% end
% 
% %
% % peakup2 =find((ddy12f(1:end-1)>=0)&(ddy12f(2:end)<0))+1; % all maxima
% % % peakdown2 =find((ddy12f(1:end-1)<=0)&(ddy12f(2:end)>0))+1;
% % peakdown2 = [];
% % 
% % peakup2_ind = find(dy12f(peakup2)>3*std(dy12f));
% % peakup2 = peakup2(peakup2_ind);
% % 
% % %%% minima
% % 
% % for st = 1:length(peakup2_ind)-1
% %     
% %     [miny, minx] = min(dy12f(peakup2(st):peakup2(st+1)));
% %     peakdown2 = [peakdown2 peakup2(st)+minx-1];    
% %     
% % end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% maxx = x12(peaksx_ind);
% maxy = y12f(peaksx_ind);
% 
% minx = x12(peakdown_ind);
% miny = y12f(peakdown_ind);
% 
% maxincx = x12(maxdy_x);
% maxincy = y12f(maxdy_x);
% 
% mindecx = x12(mindy_x);
% mindecy = y12f(mindy_x);
% 
% maxdy = dy12f(maxdy_x);
% mindy = dy12f(mindy_x);
% 
% %%%
% 
% maxx = maxx(2:end-1);
% maxy = maxy(2:end-1);
% 
% minx = minx(2:end);
% miny = miny(2:end);
% 
% maxincx = maxincx(1:end-1);
% maxincy = maxincy(1:end-1);
% 
% mindecx = mindecx(2:end);
% mindecy = mindecy(2:end);
% 
% maxdy = maxdy(1:end-1);
% mindy = mindy(2:end);
% 
% %%%
% 
% % segx = x12(maxdy_x(20):maxdy_x(23));
% % segy = y12f(maxdy_x(20):maxdy_x(23));
% 
% % l1 = length(maxincx);
% % l2 = length(mindecx);
% % 
% % if l1<l2
% %     mindecx = mindecx(1:length(maxincx));
% %     mindecy = mindecy(1:length(maxincy));
% % elseif l2<l1
% %     maxincx = maxincx(1:length(mindecx));
% %     maxincy = maxincy(1:length(mindecy));
% % end
% 
% APD = mindecx-maxincx;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% ma = peaksx_ind(2:end-1);
% mi = peakdown_ind(2:end);
% 
% % if ma(1)<mi(1)
% %     ma = ma(2:end);
% % end
% 
% repslope = [];
% repslopewin = [1 25];
% 
% for count = 1:length(ma)-1
%     
%     sl_interval = ma(count+1)-mi(count);
%         
%     b1 = round(repslopewin(1)*sl_interval/100);
%     b2 = round(repslopewin(2)*sl_interval/100);
%     
%     meansl = mean(dy12f(mi(count)+b1:mi(count)+b2));
%     repslope = [repslope meansl];
%     
% end
% 
% out.maxx = maxx;
% out.maxy = maxy;
% out.minx = minx;
% out.miny = miny;
% out.maxincx = maxincx;
% out.maxincy = maxincy;
% out.mindecx = mindecx;
% out.mindecy = mindecy;
% out.maxdy = maxdy;
% out.mindy = mindy;
% out.APD = APD;
% out.repslope = repslope;
