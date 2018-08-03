function p = poisstestu(x,y,n)
x=reshape(x,1,numel(x));
y=reshape(y,1,numel(y));

x(isnan(x))=[];
y(isnan(y))=[];
if mean(y)>=mean(x)
    diff=mean(y)-mean(x);
    dist = horzcat(x,y);
    diffshuff = zeros(n,1);
    for i = 1:n
        order = randperm(numel(dist));
        shuffle = dist(order);
        xshuff = shuffle(1:numel(x));
        yshuff = shuffle((numel(x)+1):numel(dist));
        diffshuff(i) = mean(yshuff)-mean(xshuff);
    end
else
    diff=mean(x)-mean(y);
    dist = horzcat(x,y);
    diffshuff = zeros(n,1);
    for i = 1:n
        order = randperm(numel(dist));
        shuffle = dist(order);
        xshuff = shuffle(1:numel(x));
        yshuff = shuffle((numel(x)+1):numel(dist));
        diffshuff(i) = mean(xshuff)-mean(yshuff);
    end
end
p = sum(abs(diffshuff)>=diff)/n;
end