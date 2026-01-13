clear all

%%
[v, f] = loadHull("build/hole.m");
trisurf(f, v(:,1), v(:,2), v(:,3), 'FaceColor', 'Cyan', 'FaceAlpha', 0.3, 'LineWidth', 0.001);
hold on;


P = readmatrix("build/nonmani.txt");
scatter3(v(P+1,1), v(P+1,2), v(P+1,3), "green"); hold on;



%%

[v, f] = loadHull("build/prev.m");
trisurf(f, v(:,1), v(:,2), v(:,3), 'FaceColor', 'red', 'FaceAlpha', 0.3, 'LineWidth', 0.001);
hold on;


%%
[v, f] = loadHull("build/deleted.m");
trisurf(f, v(:,1), v(:,2), v(:,3), 'FaceColor', 'red', 'FaceAlpha', 0.3, 'LineWidth', 0.001);
hold on;


h = readmatrix("build/horizon.txt");
scatter3(v(h(:),1), v(h(:),2), v(h(:),3), 'blue'); hold on;

%%
[v, f] = loadHull("build/deleted.m");
trisurf(f, v(:,1), v(:,2), v(:,3), 'FaceColor', 'red', 'FaceAlpha', 0.3, 'LineWidth', 0.001);
hold on;

h = readmatrix("build/visible.txt");
scatter3(v(h(:),1), v(h(:),2), v(h(:),3), 'blue'); hold on;



%%
a = 20;
scatter3(v(a,1), v(a,2), v(a,3), 'green'); 

%%
sum((v(21,:) - v(5,:)).^2)
%%
P = readmatrix("build/test_circum.txt");

scatter3(P(2:end,1), P(2:end,2), P(2:end,3)); hold on;
scatter3(P(1,1), P(1,2), P(1,3), "red");

%%
[vb, fb] = loadHull("build/before.m");
[va, fa] = loadHull("build/after.m");

trisurf(fb, vb(:,1), vb(:,2), vb(:,3), 'FaceColor', 'Cyan', 'FaceAlpha', 0.3, 'LineWidth', 0.001);
hold on;
trisurf(fa, va(:,1), va(:,2), va(:,3), 'FaceColor', 'red', 'FaceAlpha', 0.3, 'LineWidth', 0.001);


P = readmatrix("build/test_intersect.txt");
for i = 1:size(P, 1)
    scatter3(P(i,1), P(i,2), P(i,3));
end

%%

clear all; close all; hold on;

files = dir("build/*.m");
files = files(2:end,:);

cmap = lines(numel(files));
for k = 1:numel(files)
    [V, F] = loadHull(fullfile(files(k).folder, files(k).name));
    trisurf(F, V(:,1), V(:,2), V(:,3), ...
        'FaceAlpha', 0.5, 'FaceColor', cmap(k,:));
    %pause;
end

axis equal; view(3);



%%
[cd, fd] = loadHull("build/of2.m");
trisurf(fd, cd(:,1), cd(:,2), cd(:,3), "FaceColor", "cyan", "FaceAlpha", 0.3)







%%


deg = readmatrix("build/degenerate.txt");

scatter3(deg(1,1), deg(1,2), deg(1,3), '.', 'blue');
scatter3(deg(2,1), deg(2,2), deg(2,3), '.', 'blue');
scatter3(deg(3,1), deg(3,2), deg(3,3), '.', 'red');

%%
[cd, fd] = loadHull("build/cube.m");
trisurf(fd, cd(:,1), cd(:,2), cd(:,3), "FaceColor", "cyan", "FaceAlpha", 0.3)

%%
[cd, fd] = loadHull("build/cube_duplicates.m");
trisurf(fd, cd(:,1), cd(:,2), cd(:,3), "FaceColor", "cyan", "FaceAlpha", 0.3)

%%
[cd, fd] = loadHull("build/cube_degenerate.m");
trisurf(fd, cd(:,1), cd(:,2), cd(:,3), "FaceColor", "cyan", "FaceAlpha", 0.3)


%%

[cube1, fcube1] = loadHull("build/cube1.m");
[cube2, fcube2] = loadHull("build/cube2.m");
[I, fI] = loadHull("build/intersection_cubes.m");

trisurf(fcube1, cube1(:,1), cube1(:,2), cube1(:,3), "FaceColor", "cyan", "FaceAlpha", 0.3)

hold on
trisurf(fcube2, cube2(:,1), cube2(:,2), cube2(:,3), "FaceColor", "red", "FaceAlpha", 0.3)

trisurf(fI, I(:,1), I(:,2), I(:,3), "FaceColor", "blue", "FaceAlpha", 0.3)

%% Contained 

[cube1, fcube1] = loadHull("build/cube1.m");
[cube2, fcube2] = loadHull("build/cube3.m");
[I, fI] = loadHull("build/intersection_cubes_contained.m");

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
%trisurf(c2, L2(:,1), L2(:,2), L2(:,3), 'FaceColor', 'red', 'FaceAlpha', 0.3, 'LineWidth', 0.001);
%trisurf(c3, L3(:,1), L3(:,2), L3(:,3), 'FaceColor', 'green', 'FaceAlpha', 0.3, 'LineWidth', 0.001);


%trisurf(c4, L4(:,1), L4(:,2), L4(:,3), 'FaceColor', 'blue', 'FaceAlpha', 0.3, 'LineWidth', 0.001);
trisurf(c5, L5(:,1), L5(:,2), L5(:,3), 'FaceColor', 'yellow', 'FaceAlpha', 0.3, 'LineWidth', 0.001);

%%

%trisurf(c1, L1(:,1), L1(:,2), L1(:,3), 'FaceColor', 'Cyan', 'FaceAlpha', 0.3, 'LineWidth', 0.001);

trisurf(c2, L2(:,1), L2(:,2), L2(:,3), 'FaceColor', 'red', 'FaceAlpha', 0.3, 'LineWidth', 0.001); hold on;
trisurf(c3, L3(:,1), L3(:,2), L3(:,3), 'FaceColor', 'green', 'FaceAlpha', 0.3, 'LineWidth', 0.001);


trisurf(c4, L4(:,1), L4(:,2), L4(:,3), 'FaceColor', 'blue', 'FaceAlpha', 0.3, 'LineWidth', 0.001);
%trisurf(c5, L5(:,1), L5(:,2), L5(:,3), 'FaceColor', 'yellow', 'FaceAlpha', 0.3, 'LineWidth', 0.001);



%%
[p1, f1] = loadHull("build/problematic_pre_1.m");
[p2, f2] = loadHull("build/problematic_pre_2.m");
[p3, f3] = loadHull("build/problematic_post.m");

trisurf(f1, p1(:,1), p1(:,2), p1(:,3), 'FaceColor', 'red', 'FaceAlpha', 0.3, 'LineWidth', 0.001);
hold on;
trisurf(f2, p2(:,1), p2(:,2), p2(:,3), 'FaceColor', 'red', 'FaceAlpha', 0.3, 'LineWidth', 0.001);
trisurf(f3, p3(:,1), p3(:,2), p3(:,3), 'FaceColor', 'blue', 'FaceAlpha', 0.3, 'LineWidth', 0.001);

a = 7;
scatter3(p1(a,1), p1(a,2), p1(a,3), "black"); hold on;

%%
P = readmatrix("build/outside.txt");
scatter3(P(:,1), P(:,2), P(:,3), "green"); hold on;

P = readmatrix("build/intersect.txt");
scatter3(P(:,1), P(:,2), P(:,3), "cyan"); hold on;



%%
[v, f] = loadHull("build/of.m");
trisurf(f, v(:,1), v(:,2), v(:,3), 'FaceColor', 'red', 'FaceAlpha', 0.3, 'LineWidth', 0.001);
hold on;



%%
function [vertices, faces] = loadHull(filename)
    run(filename);
end