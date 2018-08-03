function stat = chistat_ind(tab)

rowtot=sum(tab,2);
coltot=sum(tab,1);
tot=sum(rowtot);

expect=rowtot*coltot/tot;

stat=sum(sum((tab-expect).^2./expect));

end