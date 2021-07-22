function PH = boxutil3(x,notch,lb,lf,sym,vert,whis,barva)
%BOXUTIL Produces a single box plot.
%   BOXUTIL(X) is a utility function for BOXPLOT, which calls
%   BOXUTIL once for each column of its first argument. Use
%   BOXPLOT rather than BOXUTIL. 
%   Copyright (c) 1993-96 by The MathWorks, Inc.
%   $Revision: 2.5 $  $Date: 1996/10/07 15:56:15 $
% x: data
% notch: 0-rectangular, 1-notched-box
% lb: x position
% lf: box width
% sym: symbol for the outlayer values
% vert: 0-horizontal box, 1-vertical box
% vhis

% Make sure X is a vector.
if min(size(x)) ~= 1, 
    error('Requires a vector first argument.'); 
end

if nargin ~= 8
    error('Requires eight input arguments.');
end

% define the median and the quantiles
med = prctile(x,50);
q1 = prctile(x,25);
q3 = prctile(x,75);

% find the extreme values (to determine where whiskers appear)
vhi = q3+whis*(q3-q1);
[l,ind] = min(abs((x - vhi) + (x > vhi) * max(abs(x - vhi))));
upadj = x(ind);
vlo = q1-whis*(q3-q1);
[l,ind] = min(abs((vlo - x) + (vlo > x) * max(abs(vlo - x))));
loadj = x(ind);

x1 = lb*ones(1,2);
x2 = x1+[-0.25*lf,0.25*lf];
yy = x(x<loadj | x > upadj);

if isempty(yy)
    yy = loadj;
    sym = 'g.';
end

xx = lb*ones(1,length(yy));
    lbp = lb + 0.5*lf;
    lbm = lb - 0.5*lf;

% Set up (X,Y) data for notches if desired.
if ~notch
    xx2 = [lbm lbp lbp lbm lbm];
    yy2 = [q3 q3 med med q3];
    yy3 = [med med q1 q1 med];
    lnm = lb-0.25*lf;
    lnp = lb+0.25*lf;
else
    n1 = med + 1.57*(q3-q1)/sqrt(length(x));
    n2 = med - 1.57*(q3-q1)/sqrt(length(x));
    if n1>q3, n1 = q3; end
    if n2<q1, n2 = q1; end
    lnm = lb-0.25*lf;
    lnp = lb+0.25*lf;
    xx2 = [lbm lbp lbp lnp lnm lbm lbm];
    yy2 = [q3 q3 n1 med med n1 q3];
    yy3 = [q1 q1 n2 med med n2 q1];
end

% Determine if the boxes are vertical or horizontal.
% The difference is the choice of x and y in the plot command.
if vert
    p1=patch(xx2,yy2,barva); % polnilo
    patch(xx2,yy3,barva); % polnilo
    hold on
	
	if isempty(sym)
		pout1=plot(x1,[q3 upadj],'k--',x1,[loadj q1],'k--',...
			x2,[loadj loadj],'k-',...
            x2,[upadj upadj],'k-',xx2,yy2,'k-',xx2,yy3,'k-');
    else
		pout1 = plot(x1,[q3 upadj],'k--',x1,[loadj q1],'k--',...
			x2,[loadj loadj],'k-',...
			x2,[upadj upadj],'k-',xx2,yy2,'k-',xx2,yy3,'k-',xx,yy,sym);        
    end
    set(pout1,'color','k')
    set(pout1(7),'color','w','markerfacecolor',[0.8 0.8 0.8])
	
line([lnm lnp],[med med],'color','k','linewidth',2) ;

else
    patch(yy2,xx2,barva); % polnilo
    patch(yy3,xx2,barva); % polnilo
	hold on
	
	if isempty(sym)
		pout1=plot([q3 upadj],x1,'k--',[loadj q1],x1,'k--',...
			[loadj loadj],x2,'k-',...
			[upadj upadj],x2,'k-',yy2,xx2,'k-',yy3,xx2,'k-')
	else
		pout1=plot([q3 upadj],x1,'k--',[loadj q1],x1,'k--',...
			[loadj loadj],x2,'k-',...
			[upadj upadj],x2,'k-',yy2,xx2,'k-',yy3,xx2,'k-',yy,xx,sym)
	end
	set(pout1,'color','b')
    
	line([med med],[lnm lnp],'color','k','linewidth',2) ;
	
end

PH = p1;
return


