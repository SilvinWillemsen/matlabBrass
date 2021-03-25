epsilon = 0;

if connectedWithP
    etaPrevP = (wpPrev(1) - upPrev(end));
    
    sig0 = 50;
    rForce = (1 - sig0 / k) / (1 + sig0 / k);
    oOP = (h * (1 + sig0 / k) * (1-alf)) / (2 * h * (alf + epsilon) + 2 * k^2 * (1 + sig0 / k) * (1-alf));
    
    FP = ((wpNext(1) - upNext(end)) + rForce * etaPrevP) * oOP;
    
    upNext(end) = upNext(end) + k^2/h * FP;
    wpNext(1) = wpNext(1) - k^2/h * FP;
    
    if correctV
        if correctVirtual
            etaDiv = 0.5;
            etaPrevV1 = (wvPrev(2) - uvPrev(end)) * etaDiv;
            etaPrevV2 = (wvPrev(1) - uvPrev(end-1)) * etaDiv;
        %     
            oOPV = (h * (1 + sig0 / k) * (1-alf)) / (2 * h * (alf + epsilon) + 2 * etaDiv * k^2 * (1 + sig0 / k) * (1-alf));

            FV1 = ((wvNext(2) - uvNext(end)) * etaDiv + rForce * etaPrevV1) * oOPV;
            FV2 = ((wvNext(1) - uvNext(end-1)) * etaDiv + rForce * etaPrevV2) * oOPV;

            uvNext(end) = uvNext(end) + k^2/h * FV1;
            wvNext(2) = wvNext(2) - k^2/h * FV1;

            uvNext(end-1) = uvNext(end-1) + k^2/h * FV2;
            wvNext(1) = wvNext(1) - k^2/h * FV2;
        else
            etaPrevV = (wvPrev(1) - uvPrev(end));
    
            sig0 = 50;
            rForce = (1 - sig0 / k) / (1 + sig0 / k);
            oOP = (h * (1 + sig0 / k) * alf) / (2 * h * ((1-alf) + epsilon) + 2 * k^2 * (1 + sig0 / k) * alf);

            FV = ((wvNext(1) - uvNext(end)) + rForce * etaPrevV) * oOP;

            uvNext(end) = uvNext(end) + k^2/h * FV;
            wvNext(1) = wvNext(1) - k^2/h * FV;
        end
%         if FV1 ~= 0
%             disp("wait")
%         end
    end

else
    etaPrevV = (wvPrev(1) - uvPrev(end));
    
    sig0 = 10;
    rForce = (1 - sig0 / k) / (1 + sig0 / k);
    oOP = (h * (1 + sig0 / k) * (1-alf)) / (2 * h * (alf + epsilon) + 2 * k^2 * (1 + sig0 / k) * (1-alf));
    
    FV = ((wvNext(1) - uvNext(end)) + rForce * etaPrevV) * oOP;
    
    uvNext(end) = uvNext(end) + k^2/h * FV;
    wvNext(1) = wvNext(1) - k^2/h * FV;
end
    