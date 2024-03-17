%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code creates a mac file containing a PLA slab with grid interior
% fill pattern in the Y direction (vertical printing, vertical hexagons)
% with the user's infill percentage specification.
%
% Contributors: Cláudia Reis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

slab_margin = 0.04;  %cm --> 0.02 x 2 (Simplify3D)
square_margin = 0.030; % cm --> gap --> dependes on the infill (0.025-0.035 cm)
solid_layers = 0.06; %cm --> 0.02 x 3 (top and bottom layers)

infill = input("What is the desired infill percentage?");
number_squares = (infill*13)/15; % 15% corresponds to 13 full squares horizontally
number_squares_x = (infill*2)/28; % 28% corresponds to 2 full squares vertically

if floor(number_squares)==number_squares % full number of squares
    if rem(floor(number_squares),2) == 0 % even number
        diagonal = round((((10-2*slab_margin)-(floor(number_squares)-1)*square_margin)/number_squares), 4);
        repeater_z_1 = floor(number_squares);
        repeater_z_2 = floor(number_squares) + 1; % halfs
        repeater_x_1 = floor(number_squares_x);
        repeater_x_2 = floor(number_squares_x) + 1;
    else % odd number
        diagonal = round((((10-2*slab_margin)-floor(number_squares)*square_margin)/number_squares), 4);
        repeater_z_1 = floor(number_squares) + 1; % halfs
        repeater_z_2 = floor(number_squares);
        repeater_x_1 = floor(number_squares_x) + 1;
        repeater_x_2 = floor(number_squares_x);
    end
else % not a full number of hexagons
    if rem(floor(number_squares),2) == 0 % even number
        diagonal = round((((10-2*slab_margin)-(floor(number_squares)+1)*square_margin)/number_squares), 4);
        repeater_z_1 = floor(number_squares) + 2;
        repeater_z_2 = floor(number_squares) + 1;
        repeater_x_1 = floor(number_squares_x) + 2;
        repeater_x_2 = floor(number_squares_x) + 1;
    else % odd number
        diagonal = round((((10-2*slab_margin)-floor(number_squares)*square_margin)/number_squares), 4);
        repeater_z_1 = floor(number_squares) + 1;
        repeater_z_2 = floor(number_squares) + 2;
        repeater_x_1 = floor(number_squares_x) + 1;
        repeater_x_2 = floor(number_squares_x) + 2;
    end
end

square_edge = diagonal/(sqrt(2));

%% WRITING IN MAC FILE

fid = fopen('main_grid_y.mac', 'w');

% INTRODUCTION
fprintf(fid, '\n%s \n%s %d\n', '# Grid Y direction', '# Infill percentage:', infill);

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

% square 1
fprintf(fid, '\n\n# SQUARE 1');
fprintf(fid, '\n/gate/slab/daughters/name square1');
fprintf(fid, '\n/gate/slab/daughters/insert box');
fprintf(fid, '\n%s %.2f %s', '/gate/square1/geometry/setXLength', square_edge, 'cm');
fprintf(fid, '\n%s %.2f %s', '/gate/square1/geometry/setYLength', 10-(solid_layers*2), 'cm');
fprintf(fid, '\n%s %.2f %s', '/gate/square1/geometry/setZLength', square_edge, 'cm');
fprintf(fid, '\n/gate/square1/placement/setRotationAxis 0 1 0');
fprintf(fid, '\n/gate/square1/placement/setRotationAngle 45 deg');
fprintf(fid, '\n/gate/square1/placement/setTranslation 0 0 0 cm');
fprintf(fid, '\n/gate/square1/setMaterial Air');
fprintf(fid, '\n/gate/square1/vis/setColor yellow');

% square 2
fprintf(fid, '\n\n# SQUARE 2');
fprintf(fid, '\n/gate/slab/daughters/name square2');
fprintf(fid, '\n/gate/slab/daughters/insert box');
fprintf(fid, '\n%s %.2f %s', '/gate/square2/geometry/setXLength', square_edge, 'cm');
fprintf(fid, '\n%s %.2f %s', '/gate/square2/geometry/setYLength', 10-(solid_layers*2), 'cm');
fprintf(fid, '\n%s %.2f %s', '/gate/square2/geometry/setZLength', square_edge, 'cm');
fprintf(fid, '\n/gate/square2/placement/setRotationAxis 0 1 0');
fprintf(fid, '\n/gate/square2/placement/setRotationAngle 45 deg');
fprintf(fid, '\n/gate/square2/placement/setTranslation 0 0 0 cm');
fprintf(fid, '\n/gate/square2/setMaterial Air');
fprintf(fid, '\n/gate/square2/vis/setColor white');

% repeaters 1
fprintf(fid, '\n\n# REPEATER 1');
fprintf(fid, '\n/gate/square1/repeaters/insert cubicArray');
fprintf(fid, '\n%s %.0f', '/gate/square1/cubicArray/setRepeatNumberX', repeater_x_1);
fprintf(fid, '\n/gate/square1/cubicArray/setRepeatNumberY 1');
fprintf(fid, '\n%s %.0f', '/gate/square1/cubicArray/setRepeatNumberZ', repeater_z_2);
fprintf(fid, '\n%s %.3f %s %.3f %s', '/gate/square1/cubicArray/setRepeatVector', diagonal+square_margin, '0', diagonal+square_margin, 'cm');

% repeaters 2
fprintf(fid, '\n\n# REPEATER 2');
fprintf(fid, '\n/gate/square2/repeaters/insert cubicArray');
fprintf(fid, '\n%s %.0f', '/gate/square2/cubicArray/setRepeatNumberX', repeater_x_2);
fprintf(fid, '\n/gate/square2/cubicArray/setRepeatNumberY 1');
fprintf(fid, '\n%s %.0f', '/gate/square2/cubicArray/setRepeatNumberZ', repeater_z_1);
fprintf(fid, '\n%s %.3f %s %.3f %s', '/gate/square2/cubicArray/setRepeatVector', diagonal+square_margin, '0', diagonal+square_margin, 'cm');

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
fprintf(fid, '\n/gate/physics/Gamma/SetCutInRegion square1 0.01 mm');
fprintf(fid, '\n/gate/physics/Electron/SetCutInRegion square1 0.01 mm');
fprintf(fid, '\n/gate/physics/Positron/SetCutInRegion square1 0.01 mm');
fprintf(fid, '\n/gate/physics/Gamma/SetCutInRegion square2 0.01 mm');
fprintf(fid, '\n/gate/physics/Electron/SetCutInRegion square2 0.01 mm');
fprintf(fid, '\n/gate/physics/Positron/SetCutInRegion square2 0.01 mm');
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
fprintf(fid, '\n/gate/actor/stat/save output/stats_grid_y_%0.0f.txt', infill);

% EmCalculator
fprintf(fid, '\n\n#=====================================================\n# EmCalculator actor \n#=====================================================\n');
fprintf(fid, '\n/gate/actor/addActor EmCalculatorActor EM-properties');
fprintf(fid, '\n/gate/actor/EM-properties/setParticleName proton');
fprintf(fid, '\n/gate/actor/EM-properties/setEnergy 210 MeV');
fprintf(fid, '\n/gate/actor/EM-properties/save output/EMprops_grid_y_%0.0f.txt', infill);

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
fprintf(fid, '\n/gate/actor/doseDistribution/save output/output_grid_y_%0.0f.mhd', infill);

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