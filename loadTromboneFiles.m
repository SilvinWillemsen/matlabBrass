myPath = "/Users/SilvinW/repositories/brass/cppBrass/Trombone/Builds/MacOSX/build/" + mode;

if ~onlyLoadOutput
    pStateString = myPath + "/pState.csv";
    pState = load(pStateString);
    disp("pState Loaded")

    vStateString = myPath + "/vState.csv";
    vState = load(vStateString);
    disp("vState Loaded")

    alfSaveString = myPath + "/alfSave.csv";
    alfSave = load(alfSaveString);
    disp("alfSave Loaded")

    massStateString = myPath + "/massState.csv";
    massState = load(massStateString);
    disp("massState Loaded")

    MString = myPath + "/MSave.csv";
    M = load(MString);
    disp("M Loaded")

    MwString = myPath + "/MwSave.csv";
    Mw = load(MwString);
    disp("Mw Loaded")

    energyString = myPath + "/energySave.csv";
    energy = load(energyString);
    disp("energy Loaded")

    scaledTotEnergyString = myPath + "/scaledTotEnergySave.csv";
    scaledTotEnergy = load(scaledTotEnergyString);
    disp("scaledTotEnergy Loaded")

    maxMString = myPath + "/maxM.csv";
    maxM = load(maxMString);
    disp("maxM Loaded")

    maxMwString = myPath + "/maxMw.csv";
    maxMw = load(maxMwString);
    disp("maxMw Loaded")

    SSaveString = myPath + "/SSave.csv";
    Ssave = load(SSaveString);
    disp("Ssave Loaded")   
    
    statesSaveString = myPath + "/statesSave.csv";
    statesSaveCPP = load(statesSaveString);
    disp("statesSaveCPP Loaded")

end
outputString = myPath + "/output.csv";
output = load(outputString);
disp("output Loaded")

