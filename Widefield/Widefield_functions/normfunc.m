function x = normfunc(x)
    x = x - min(x(:));
    x = x./max(x(:));
end