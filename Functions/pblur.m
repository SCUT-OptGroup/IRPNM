function y = pblur(x,p,q)
if isvector(x)
    x = reshape(x,[p q]);
end
y = imgaussfilt(x,4,'FilterSize',9,'Padding','symmetric');
y = y(:);
end
