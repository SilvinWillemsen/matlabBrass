epsilon = 0;
sig0 = 10;

if correctV == 1
    etaDiv = 0.5;
    etaPrevV1 = (wvPrev(2) - uvPrev(end)) * etaDiv;
    etaPrevV2 = (wvPrev(1) - uvPrev(end-1)) * etaDiv;
    rForce = (1 - sig0 / k) / (1 + sig0 / k);
    oOPV = (h * (1 + sig0 / k) * (1-alf)) / (2 * h * (alf + epsilon) + 2 * etaDiv * k^2 * (1 + sig0 / k) * (1-alf));

    FV1 = ((wvNext(2) - uvNext(end)) * etaDiv + rForce * etaPrevV1) * oOPV;
    FV2 = ((wvNext(1) - uvNext(end-1)) * etaDiv + rForce * etaPrevV2) * oOPV;

    uvNext(end) = uvNext(end) + k^2/h * FV1;
    wvNext(2) = wvNext(2) - k^2/h * FV1;

    uvNext(end-1) = uvNext(end-1) + k^2/h * FV2;
    wvNext(1) = wvNext(1) - k^2/h * FV2;
elseif correctV == 2
    etaPrevV = (wvPrev(1) - uvPrev(end));

    rForce = (1 - sig0 / k) / (1 + sig0 / k);
    oOPV = (h * (1 + sig0 / k) * alf) / (2 * h * ((1-alf) + epsilon) + 2 * k^2 * (1 + sig0 / k) * alf);

    FV = ((wvNext(1) - uvNext(end)) + rForce * etaPrevV) * oOPV;

    uvNext(end) = uvNext(end) + k^2/h * FV;
    wvNext(1) = wvNext(1) - k^2/h * FV;
elseif correctV == 3
    Iu = zeros(1, length(uv));
    Iu(end-1:end) = [0.5, 0.5];
    Ju = 1/h * Iu';

    Iw = zeros(1, length(wv));
    Iw(1:2) = [0.5, 0.5];
    Jw = 1/h * Iw';

    etaPrevV = Iw * wvPrev - Iu * uvPrev;
    beta = 1000000;
    rForce = (1 - sig0 / k) / (1 + sig0 / k);
    oOPV = ((1 + sig0 / k) * beta) / (2 + k^2 * (1 + sig0 / k) * beta * (Iw * Jw + Iu * Ju));

    FV = ((Iw * wvNext - Iu * uvNext) + rForce * etaPrevV) * oOPV;

    uvNext = uvNext + Ju * k^2 * FV;
    wvNext = wvNext - Jw * k^2 * FV;

end