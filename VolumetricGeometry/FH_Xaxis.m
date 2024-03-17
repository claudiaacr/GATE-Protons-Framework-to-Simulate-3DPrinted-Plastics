%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code creates a mac file containing a PLA slab with full honeycomb
% interior fill pattern in the X direction (horizontal printing) with the
% user's infill percentage specification.
%
% Contributors: Cláudia Reis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

slab_margin = 0.04;  %cm --> 0.02 x 2 (Simplify3D)
hex_margin = 0.04;  %cm --> gap (may change)
solid_layers = 0.06; %cm --> 0.02 x 3 (top and bottom layers)

infill = input('What is the desired infill percentage? ');
% 13% infill corresponds to 12 full hexagons vertially
number_hexs = (infill*12)/13;

if rem(floor(number_hexs),2) == 0 % even number
    radius = round((((10-floor(number_hexs)*hex_margin)/number_hexs)/2), 4);
else % odd number
    radius = round((((10-(floor(number_hexs)+1)*hex_margin)/number_hexs)/2), 4);
end

dist_corners = round(1.154701 * radius * 2, 2); %cm 
hex_edge = round(2 * (radius/tan(pi/3)), 2);

%% WRITING THE MAC FILE

fid = fopen('main_hexs_x.mac', 'w');

% INTRODUCTION
fprintf(fid, '\n%s \n%s %d\n\n', '# Full Honeycomb X direction', '# Infill percentage:', infill);

% VERBOSE AND VISUALISATION
fprintf(fid, '#=====================================================\n# VERBOSE and VISUALISATION\n#=====================================================\n\n/control/execute mac/verbose.mac\n/control/execute mac/visu.mac\n');

% GEOMETRY
fprintf(fid, '\n#=====================================================\n# GEOMETRY\n#=====================================================\n');
fprintf(fid, '\n/gate/geometry/setMaterialDatabase data/GateMaterials.db\n/gate/geometry/setIonisationPotential Water 78 eV\n/gate/geometry/setIonisationPotential PLA100 77.958580668515992 eV');

% world
fprintf(fid, '\n\n# WORLD');
fprintf(fid, '\n/gate/world/setMaterial Air');
fprintf(fid, '\n/gate/world/geometry/setXLength 5.0 m');
fprintf(fid, '\n/gate/world/geometry/setYLength 5.0 m');
fprintf(fid, '\n/gate/world/geometry/setZLength 5.0 m');

% slab
fprintf(fid, '\n\n# SLAB');
fprintf(fid, '\n/gate/world/daughters/name slab');
fprintf(fid, '\n/gate/world/daughters/insert box');
fprintf(fid, '\n/gate/slab/setMaterial PLA100');
fprintf(fid, '\n/gate/slab/geometry/setXLength 1 cm');
fprintf(fid, '\n/gate/slab/geometry/setYLength 10 cm');
fprintf(fid, '\n/gate/slab/geometry/setZLength 10 cm');
fprintf(fid, '\n/gate/slab/placement/setTranslation 0 0 0 cm');
fprintf(fid, '\n/gate/slab/vis/setColor magenta');

% hexagon 1
fprintf(fid, '\n\n# HEXAGON 1');
fprintf(fid, '\n/gate/slab/daughters/name hexagon1');
fprintf(fid, '\n/gate/slab/daughters/insert hexagone');
fprintf(fid, '\n%s %.4f %s', '/gate/hexagon1/geometry/setRadius', radius, 'cm');
fprintf(fid, '\n%s %.4f %s', '/gate/hexagon1/geometry/setHeight', 1-(solid_layers*2), 'cm');
fprintf(fid, '\n/gate/hexagon1/placement/setRotationAxis 0 1 0');
fprintf(fid, '\n/gate/hexagon1/placement/setRotationAngle 90 deg');
fprintf(fid, '\n/gate/hexagon1/placement/setTranslation 0 0 0 cm');
fprintf(fid, '\n/gate/hexagon1/setMaterial Air');
fprintf(fid, '\n/gate/hexagon1/vis/setColor yellow');

% hexagon 2
fprintf(fid, '\n\n# HEXAGON 2');
fprintf(fid, '\n/gate/slab/daughters/name hexagon2');
fprintf(fid, '\n/gate/slab/daughters/insert hexagone');
fprintf(fid, '\n%s %.4f %s', '/gate/hexagon2/geometry/setRadius', radius, 'cm');
fprintf(fid, '\n%s %.4f %s', '/gate/hexagon2/geometry/setHeight', 1-(solid_layers*2), 'cm');
fprintf(fid, '\n/gate/hexagon2/placement/setRotationAxis 0 1 0');
fprintf(fid, '\n/gate/hexagon2/placement/setRotationAngle 90 deg');
fprintf(fid, '\n/gate/hexagon2/placement/setTranslation 0 0 0 cm');
fprintf(fid, '\n/gate/hexagon2/setMaterial Air');
fprintf(fid, '\n/gate/hexagon2/vis/setColor yellow');

% repeaters 1
fprintf(fid, '\n\n# REPEATER 1');
fprintf(fid, '\n/gate/hexagon1/repeaters/insert cubicArray');
fprintf(fid, '\n/gate/hexagon1/cubicArray/setRepeatNumberX 1');

if floor(number_hexs)==number_hexs % full number of hexagons
    if rem(floor(number_hexs),2) == 0 % even number
        repeater_y_1 = number_hexs + 1; % two halfs
        repeater_y_2 = repeater_y_1 - 1;
    else % odd number
        repeater_y_1 = number_hexs; % complete hexagons
        repeater_y_2 = repeater_y_1 + 1;
    end
else % not a full number of hexagons
    if rem(floor(number_hexs),2) == 0 % even number
        repeater_y_1 = floor(number_hexs)+1;
        repeater_y_2 = repeater_y_1 + 1;
    else % odd number
        repeater_y_1 = floor(number_hexs)+2;
        repeater_y_2 = repeater_y_1 - 1;
    end
end

if infill > 35
    fprintf(fid, '\n%s %.0f', '/gate/hexagon1/cubicArray/setRepeatNumberY', repeater_y_1-2);
else
    fprintf(fid, '\n%s %.0f', '/gate/hexagon1/cubicArray/setRepeatNumberY', repeater_y_1);
end

half_slab_edge = 5 - (dist_corners/2) - hex_margin;
count_hex_edge = 0;
count_dist_corners = 0;
a = 0;

while half_slab_edge >= dist_corners && a < 100 % bigger than hex_edge
    if half_slab_edge >= (hex_edge + hex_margin) % the next one after the central hexagon
        half_slab_edge = half_slab_edge - hex_edge - hex_margin;
        count_hex_edge = count_hex_edge + 1;
    end
    if half_slab_edge >= (dist_corners + hex_margin)
        half_slab_edge = half_slab_edge - dist_corners - hex_margin;
        count_dist_corners = count_dist_corners + 1;
    end
    a = a + 1;
end
    
if count_hex_edge == count_dist_corners % then comes hex_edge or nothing (complete hexagon) --> doesn't start a new hexagon
    repeater_z_1 = 2 * count_dist_corners + 1; % each side + central hexagon
    repeater_z_2 = 2 * (count_hex_edge + 1);
else % dist_corners comes (incomplete, but counts) or it doesn't 
    if half_slab_edge == 0 % complete
        repeater_z_1 = 2 * count_dist_corners + 1;
    else % incomplete
        repeater_z_1 = 2 * (count_dist_corners + 1) + 1;
    repeater_z_2 = 2 * count_hex_edge;
    end
end

fprintf(fid, '\n%s %.0f', '/gate/hexagon1/cubicArray/setRepeatNumberZ', repeater_z_1);
fprintf(fid, '\n%s %.5f %.5f %s', '/gate/hexagon1/cubicArray/setRepeatVector 0', radius*2+hex_margin, dist_corners+hex_edge+2*+hex_margin, 'cm');

% repeaters 2
fprintf(fid, '\n\n# REPEATER 2');
fprintf(fid, '\n/gate/hexagon2/repeaters/insert cubicArray');
fprintf(fid, '\n/gate/hexagon2/cubicArray/setRepeatNumberX 1');

if infill > 35
    fprintf(fid, '\n%s %.0f', '/gate/hexagon2/cubicArray/setRepeatNumberY', repeater_y_2-2);
else
    fprintf(fid, '\n%s %.0f', '/gate/hexagon2/cubicArray/setRepeatNumberY', repeater_y_2);
end

fprintf(fid, '\n%s %.0f', '/gate/hexagon2/cubicArray/setRepeatNumberZ', repeater_z_2);
fprintf(fid, '\n%s %.5f %.5f %s', '/gate/hexagon2/cubicArray/setRepeatVector 0', radius*2+hex_margin, dist_corners+hex_edge+2*+hex_margin, 'cm');

% margin 1 (up)
fprintf(fid, '\n\n# MARGIN 1');
fprintf(fid, '\n/gate/slab/daughters/name margin1');
fprintf(fid, '\n/gate/slab/daughters/insert box');
fprintf(fid, '\n/gate/margin1/setMaterial PLA100');
fprintf(fid, '\n/gate/margin1/geometry/setXLength 1 cm');
fprintf(fid, '\n%s %.2f %s', '/gate/margin1/geometry/setYLength', slab_margin, 'cm');
fprintf(fid, '\n/gate/margin1/geometry/setZLength 10 cm');
fprintf(fid, '\n%s %.2f %s', '/gate/margin1/placement/setTranslation 0', 5-(slab_margin/2), '0 cm');
fprintf(fid, '\n/gate/margin1/vis/setColor white');

% margin 2 (bottom)
fprintf(fid, '\n\n# MARGIN 2');
fprintf(fid, '\n/gate/slab/daughters/name margin2');
fprintf(fid, '\n/gate/slab/daughters/insert box');
fprintf(fid, '\n/gate/margin2/setMaterial PLA100');
fprintf(fid, '\n/gate/margin2/geometry/setXLength 1 cm');
fprintf(fid, '\n%s %.2f %s', '/gate/margin2/geometry/setYLength', slab_margin, 'cm');
fprintf(fid, '\n/gate/margin2/geometry/setZLength 10 cm');
fprintf(fid, '\n%s %.2f %s', '/gate/margin2/placement/setTranslation 0', -(5-(slab_margin/2)), '0 cm');
fprintf(fid, '\n/gate/margin2/vis/setColor white');

% margin 3 (left)
fprintf(fid, '\n\n# MARGIN 3');
fprintf(fid, '\n/gate/slab/daughters/name margin3');
fprintf(fid, '\n/gate/slab/daughters/insert box');
fprintf(fid, '\n/gate/margin3/setMaterial PLA100');
fprintf(fid, '\n/gate/margin3/geometry/setXLength 1 cm');
fprintf(fid, '\n/gate/margin3/geometry/setYLength 10 cm');
fprintf(fid, '\n%s %.2f %s', '/gate/margin3/geometry/setZLength', slab_margin, 'cm');
fprintf(fid, '\n%s %.2f %s', '/gate/margin3/placement/setTranslation 0 0', -(5-(slab_margin/2)), 'cm');
fprintf(fid, '\n/gate/margin3/vis/setColor white');

% margin 4 (right)
fprintf(fid, '\n\n# MARGIN 4');
fprintf(fid, '\n/gate/slab/daughters/name margin4');
fprintf(fid, '\n/gate/slab/daughters/insert box');
fprintf(fid, '\n/gate/margin4/setMaterial PLA100');
fprintf(fid, '\n/gate/margin4/geometry/setXLength 1 cm');
fprintf(fid, '\n/gate/margin4/geometry/setYLength 10 cm');
fprintf(fid, '\n%s %.2f %s', '/gate/margin4/geometry/setZLength', slab_margin, 'cm');
fprintf(fid, '\n%s %.2f %s', '/gate/margin4/placement/setTranslation 0 0', 5-(slab_margin/2), 'cm');
fprintf(fid, '\n/gate/margin4/vis/setColor white');

% PHYSICS
fprintf(fid, '\n\n#=====================================================\n# PHYSICS \n#=====================================================\n');
fprintf(fid, '\n/gate/physics/addPhysicsList QGSP_BIC_EMZ');
fprintf(fid, '\n/gate/physics/Gamma/SetCutInRegion world 1 mm');
fprintf(fid, '\n/gate/physics/Electron/SetCutInRegion world 1 mm');
fprintf(fid, '\n/gate/physics/Positron/SetCutInRegion world 1 mm');
fprintf(fid, '\n/gate/physics/Gamma/SetCutInRegion slab 0.01 mm');
fprintf(fid, '\n/gate/physics/Electron/SetCutInRegion slab 0.01 mm');
fprintf(fid, '\n/gate/physics/Positron/SetCutInRegion slab 0.01 mm');
fprintf(fid, '\n/gate/physics/Gamma/SetCutInRegion hexagon1 0.01 mm');
fprintf(fid, '\n/gate/physics/Electron/SetCutInRegion hexagon1 0.01 mm');
fprintf(fid, '\n/gate/physics/Positron/SetCutInRegion hexagon1 0.01 mm');
fprintf(fid, '\n/gate/physics/Gamma/SetCutInRegion hexagon2 0.01 mm');
fprintf(fid, '\n/gate/physics/Electron/SetCutInRegion hexagon2 0.01 mm');
fprintf(fid, '\n/gate/physics/Positron/SetCutInRegion hexagon2 0.01 mm');
fprintf(fid, '\n/gate/physics/Gamma/SetCutInRegion margin1 0.01 mm');
fprintf(fid, '\n/gate/physics/Electron/SetCutInRegion margin1 0.01 mm');
fprintf(fid, '\n/gate/physics/Positron/SetCutInRegion margin1 0.01 mm');
fprintf(fid, '\n/gate/physics/Gamma/SetCutInRegion margin2 0.01 mm');
fprintf(fid, '\n/gate/physics/Electron/SetCutInRegion margin2 0.01 mm');
fprintf(fid, '\n/gate/physics/Positron/SetCutInRegion margin2 0.01 mm');
fprintf(fid, '\n/gate/physics/Gamma/SetCutInRegion margin3 0.01 mm');
fprintf(fid, '\n/gate/physics/Electron/SetCutInRegion margin3 0.01 mm');
fprintf(fid, '\n/gate/physics/Positron/SetCutInRegion margin3 0.01 mm');
fprintf(fid, '\n/gate/physics/Gamma/SetCutInRegion margin4 0.01 mm');
fprintf(fid, '\n/gate/physics/Electron/SetCutInRegion margin4 0.01 mm');
fprintf(fid, '\n/gate/physics/Positron/SetCutInRegion margin4 0.01 mm');
fprintf(fid, '\n\n/gate/physics/displayCuts');

% STATISTICS
fprintf(fid, '\n\n#=====================================================\n# Statistics actor \n#=====================================================\n');
fprintf(fid, '\n/gate/actor/addActor SimulationStatisticActor stat');
fprintf(fid, '\n/gate/actor/stat/saveEveryNSeconds 120');
fprintf(fid, '\n/gate/actor/stat/save output/stats_hexs_x_%0.0f.txt', infill);

% EmCalculator
fprintf(fid, '\n\n#=====================================================\n# EmCalculator actor \n#=====================================================\n');
fprintf(fid, '\n/gate/actor/addActor EmCalculatorActor EM-properties');
fprintf(fid, '\n/gate/actor/EM-properties/setParticleName proton');
fprintf(fid, '\n/gate/actor/EM-properties/setEnergy 210 MeV');
fprintf(fid, '\n/gate/actor/EM-properties/save output/EMprops_hexs_x_%0.0f.txt', infill);

% DOSE
fprintf(fid, '\n\n#=====================================================\n# Dose actor \n#=====================================================\n');
fprintf(fid, '\n/gate/actor/addActor DoseActor doseDistribution');
fprintf(fid, '\n/gate/actor/doseDistribution/attachTo slab'); % change
fprintf(fid, '\n/gate/actor/doseDistribution/stepHitType random');
fprintf(fid, '\n/gate/actor/doseDistribution/setVoxelSize 1 1 1 mm');
fprintf(fid, '\n/gate/actor/doseDistribution/saveEveryNSeconds 30');
fprintf(fid, '\n/gate/actor/doseDistribution/enableUncertaintyDose true');
fprintf(fid, '\n/gate/actor/doseDistribution/enableEdep true');
fprintf(fid, '\n/gate/actor/doseDistribution/enableDose true');
fprintf(fid, '\n/gate/actor/doseDistribution/enableSquaredDose false');
fprintf(fid, '\n/gate/actor/doseDistribution/enableDoseToWater false');
fprintf(fid, '\n/gate/actor/doseDistribution/save output/output_hexs_x_%0.0f.mhd', infill);

% INITIALISATION
fprintf(fid, '\n\n#=====================================================\n# INITIALISATION \n#=====================================================\n');
fprintf(fid, '\n/gate/run/initialize');

% SOURCE
fprintf(fid, '\n\n#=====================================================\n# SOURCE \n#=====================================================\n');
fprintf(fid, '\n/gate/source/addSource PBS TPSPencilBeam');
fprintf(fid, '\n/gate/source/PBS/setTestFlag false');
fprintf(fid, '\n/gate/source/PBS/setPlan data/plan_4x4.txt'); % or plan_singlePB.txt
fprintf(fid, '\n/gate/source/PBS/setFlatGenerationFlag false');
fprintf(fid, '\n/gate/source/PBS/setSpotIntensityAsNbIons false');
fprintf(fid, '\n/gate/source/PBS/setSortedSpotGenerationFlag false');
fprintf(fid, '\n/gate/source/PBS/setSigmaEnergyInMeVFlag false');
fprintf(fid, '\n/gate/source/PBS/setSourceDescriptionFile data/uclhospital_bmv3.txt');

% START BEAMS
fprintf(fid, '\n\n#=====================================================\n# START BEAMS \n#=====================================================\n');
fprintf(fid, '\n/gate/random/setEngineName MersenneTwister');
fprintf(fid, '\n/gate/random/setEngineSeed auto');
fprintf(fid, '\n\n/gate/application/setTotalNumberOfPrimaries 1000000');
fprintf(fid, '\n/gate/application/start');

fclose(fid);