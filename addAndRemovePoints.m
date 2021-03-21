%% Add and remove points
if abs(N - NPrev) > 1
    disp('too fast')
end

useAvg = true;

% add point if N^n > N^{n-1}
if N > NPrev
    customIp = [alf * (alf + 1) / -((alf + 2) * (alf + 3)); ...
        2 * alf / (alf + 2); ...
        2 / (alf + 2); ...
        2 * alf / -((alf + 3) * (alf + 2))]';
    
    if ~connectedWithP % if the new N connects at v prepare the statevectors (by adding to v)
        upMp1Prev = customIp * [upPrev(end-1:end); wpPrev(1:2)];
        wpm1Prev = fliplr(customIp) * [upPrev(end-1:end); wpPrev(1:2)];
      
        if useAvg
            vNextAvg = (uvNextMph + wvNextmh) * 0.5;
            vAvg = (uvMph + wvmh) * 0.5;

            uvNext = [uvNext; vNextAvg];
            uv = [uv; vAvg];

            wvNext = [vNextAvg; wvNext];
            wv = [vAvg; wv];
        else
            uvNext = [uvNext; uvNextMph];
            uv = [uv; uvMph];

            wvNext = [wvNextmh; wvNext];
            wv = [wvmh; wv];
        end
        justShiftedToConnectedV = true;
    else % if the new N connects at p prepare the statevectors (by adding to p)
        
        %% Here, also calculate p^n as that hasn't been done yet when switching to this setting
        uvMph = uv(end) * quadIp(3) + wv(1) * quadIp(2) + wv(2) * quadIp(1);
        wvmh = uv(end-1) * quadIp(1) + uv(end) * quadIp(2) + wv(1) * quadIp(3);

        %% Calculate p^n
        up(upRange) = upPrev(upRange) - rho * c * lambda ./ SBar(upRange) .* (SHalf(upRange) .* uv(upRange) - SHalf(upRange-1) .* uv(upRange-1));
        upMp1 = upMp1Prev - rho * c * lambda / SBar(length(up)) * (SHalf(length(up)) .* uvMph - SHalf(length(up) - 1) .* uv(end));
        up(1) = upPrev(1) - rho * c * lambda / SBar(1) .* (-2 * (Ub + Ur) + 2 * SHalf(1) * uv(1));

        wp(wpRange-1) = wpPrev(wpRange-1) - rho * c * lambda ./ SBar(wpRange + length(up) - 2) .* (SHalf(wpRange + length(up) - 2) .* wv(wpRange) - SHalf(wpRange + length(up) - 3) .* wv(wpRange-1));
        wpm1 = wpm1Prev - rho * c * lambda ./ SBar(length(up)) .* (SHalf(length(up)) .* wv(1) - SHalf(length(up) - 1) .* wvmh);
        
        if radiation
            wp(end) = ((1 - rho * c * lambda * z3) * wpPrev(end) - 2 * rho * c * lambda * (v1 + z4 * p1 - (SHalf(end) .* wv(end))/SBar(end))) / (1 + rho * c * lambda * z3);
        end
        
        v1Next = v1 + k / (2 * Lr) * (wp(end) + wpPrev(end));
        p1Next = z1 / 2 * (wp(end) + wpPrev(end)) + z2 * p1;
        
        v1 = v1Next;
        p1 = p1Next;
        
        uvMph = customIp * [uv(end-1:end); wv(1:2)];
        wvmh = fliplr(customIp) * [uv(end-1:end); wv(1:2)];
        
        if useAvg
            pAvg = (upMp1 + wpm1) * 0.5;
            pPrevAvg = (upMp1 + wpm1) * 0.5;
            
            upNext = [upNext; pAvg];
            up = [up; pAvg];
            upPrev = [upPrev; pPrevAvg];
        
            wpNext = [pAvg; wpNext];
            wp = [pAvg; wp];
            wpPrev = [pPrevAvg; wpPrev];
        else
            upNext = [upNext; upMp1];
            up = [up; upMp1];
            upPrev = [upPrev; upMp1Prev];

            wpNext = [wpm1; wpNext];
            wp = [wpm1; wp];
            wpPrev = [wpm1Prev; wpPrev];
        end
    end
    [S, SHalf, SBar] = setTube(N+1, NnonExtended, n, setToOnes);
    
    disp("point added " + num2str(alf) + " " + num2str(N))
%     plotStatesAfterAddingPoint;
end

% remove point if N^n < N^{n-1}
if N < NPrev
    if connectedWithP % remove from v to be connected at p
        uvNext = uvNext(1:end-1);
        uv = uv(1:end-1);
        wvNext = wvNext(2:end);
        wv = wv(2:end);
    else  % remove from p to be connected at v
        upNext = upNext(1:end-1);
        up = up(1:end-1);
        upPrev = upPrev(1:end-1);
        
        wpNext = wpNext(2:end);
        wp = wp(2:end);
        wpPrev = wpPrev(2:end);
        
    end
    if flag
        disp("point removed")
    end
    [S, SHalf, SBar] = setTube(N+1, NnonExtended, n);
        
end
if connectedWithP
    upRange = 2:length(up)-1;         % range without boundaries
    wpRange = 2:length(wp)-1;
else
    upRange = 2:length(up);         % range without boundaries
    wpRange = 2:length(wp);
end


