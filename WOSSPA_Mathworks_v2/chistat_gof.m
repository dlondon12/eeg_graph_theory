function stat = chistat_gof(s,n,p_exp)

obs=[s,n-s];
expect=[p_exp*n,(1-p_exp)*n];

stat=sum((obs-expect).^2./expect);

end