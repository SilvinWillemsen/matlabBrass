epsilon = 0;
sig0 = 10;

if connectedWithP
    etaPrevP = (wpPrev(1) - upPrev(end));
    
    rForce = (1 - sig0 / k) / (1 + sig0 / k);
    oOP = (h * (1 + sig0 / k) * (1-alf)) / (2 * h * (alf + epsilon) + 2 * k^2 * (1 + sig0 / k) * (1-alf));
    
    FP = ((wpNext(1) - upNext(end)) + rForce * etaPrevP) * oOP;
    
    upNext(end) = upNext(end) + k^2/h * FP;
    wpNext(1) = wpNext(1) - k^2/h * FP;
else
    etaPrevV = (wvPrev(1) - uvPrev(end));
    
    rForce = (1 - sig0 / k) / (1 + sig0 / k);
    oOP = (h * (1 + sig0 / k) * (1-alf)) / (2 * h * (alf + epsilon) + 2 * k^2 * (1 + sig0 / k) * (1-alf));
    
    FV = ((wvNext(1) - uvNext(end)) + rForce * etaPrevV) * oOP;
    
    uvNext(end) = uvNext(end) + k^2/h * FV;
    wvNext(1) = wvNext(1) - k^2/h * FV;
end
    