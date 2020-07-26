function [] = myboxplot(medianvec,q1vec,q3vec,minvec,maxvec)
% MYBOXPLOT: creates a box-plot 
nd = length(medianvec);
hm  = plot(medianvec,'k*');
set(hm,'MarkerSize',15)
hold on
d = 0.4;
for k=1:nd
    % plot box
    box_x = [k-d/2 k+d/2 k+d/2 k-d/2 k-d/2];
    box_y = [q1vec(k) q1vec(k) q3vec(k) q3vec(k) q1vec(k)];
    hb = plot(box_x, box_y,'k');
    % plot whiskers
    iqr = q3vec(k)- q1vec(k);
    uw_x = [k k k-d/2 k+d/2];
    yval = min(q3vec(k)+3/2*iqr,maxvec(k));
    uw_y = [q3vec(k), yval, yval, yval];
    hw1 = plot(uw_x,uw_y);
    lw_x = [k k k-d/2 k+d/2];
    yval = max(q1vec(k)-3/2*iqr,minvec(k));
    lw_y = [q1vec(k), yval, yval, yval];
    hw2 = plot(lw_x,lw_y);
end
end

