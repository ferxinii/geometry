clear all

[cube1, fcube1] = loadHull("build/cube1.m");
[cube2, fcube2] = loadHull("build/cube2.m");
[I, fI] = loadHull("build/intersection_cubes.m");

trisurf(fcube1, cube1(:,1), cube1(:,2), cube1(:,3), "FaceColor", "cyan", "FaceAlpha", 0.3)

hold on
trisurf(fcube2, cube2(:,1), cube2(:,2), cube2(:,3), "FaceColor", "red", "FaceAlpha", 0.3)

trisurf(fI, I(:,1), I(:,2), I(:,3), "FaceColor", "blue", "FaceAlpha", 0.3)

%% Lobes 2, 3
[L2, fL2] = loadHull("build/L2.m");
[L3, fL3] = loadHull("build/L3.m");
[I, fI] = loadHull("build/I_L2_L3.m");

trisurf(fL2, L2(:,1), L2(:,2), L2(:,3), "FaceColor", "cyan", "FaceAlpha", 0.3)
hold on
trisurf(fL3, L3(:,1), L3(:,2), L3(:,3), "FaceColor", "red", "FaceAlpha", 0.3)
trisurf(fI, I(:,1), I(:,2), I(:,3), "FaceColor", "blue", "FaceAlpha", 0.3)

%% Removed intersection
[L2p, fL2p] = loadHull("build/L2_post.m");
[L3p, fL3p] = loadHull("build/L3_post.m");

%trisurf(fL2, L2(:,1), L2(:,2), L2(:,3), "FaceColor", "green", "FaceAlpha", 0.3)
%hold on
trisurf(fL2p, L2p(:,1), L2p(:,2), L2p(:,3), "FaceColor", "cyan", "FaceAlpha", 0.3)
hold on
trisurf(fL3p, L3p(:,1), L3p(:,2), L3p(:,3), "FaceColor", "red", "FaceAlpha", 0.3)

%trisurf(fL3, L3(:,1), L3(:,2), L3(:,3), "FaceColor", "green", "FaceAlpha", 0.3)

%%
L1 = readmatrix("lobes/L1.txt");
L2 = readmatrix("lobes/L2.txt");
L3 = readmatrix("lobes/L3.txt");
L4 = readmatrix("lobes/L4.txt");
L5 = readmatrix("lobes/L5.txt");

[c1, a1] = convhull(L1);
[c2, a2] = convhull(L2);
[c3, a3] = convhull(L3);
[c4, a4] = convhull(L4);
[c5, a5] = convhull(L5);

sum_vols = a1 + a2 + a3 + a4 + a5

trisurf(c1, L1(:,1), L1(:,2), L1(:,3), 'FaceColor', 'Cyan', 'FaceAlpha', 0.3, 'LineWidth', 0.001);
hold on;
trisurf(c2, L2(:,1), L2(:,2), L2(:,3), 'FaceColor', 'red', 'FaceAlpha', 0.3, 'LineWidth', 0.001);
trisurf(c3, L3(:,1), L3(:,2), L3(:,3), 'FaceColor', 'green', 'FaceAlpha', 0.3, 'LineWidth', 0.001);
trisurf(c4, L4(:,1), L4(:,2), L4(:,3), 'FaceColor', 'blue', 'FaceAlpha', 0.3, 'LineWidth', 0.001);
trisurf(c5, L5(:,1), L5(:,2), L5(:,3), 'FaceColor', 'yellow', 'FaceAlpha', 0.3, 'LineWidth', 0.001);



%%
function [vertices, faces] = loadHull(filename)
    run(filename);
end