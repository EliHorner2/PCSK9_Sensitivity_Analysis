function KLD = KLdiv(ObservedOuts, FullModelOuts)
    evalPts = linspace(min(FullModelOuts), max(FullModelOuts), 500);
    [ft, xt] = ksdensity(ObservedOuts, evalPts);
    [fc, xc] = ksdensity(FullModelOuts, evalPts);

    KLD = 0;
    for i = 1:length(ft)
        KLD = KLD + (fc(i)*log(fc(i)/ft(i)))*(max(FullModelOuts) - min(FullModelOuts)/500);
    end