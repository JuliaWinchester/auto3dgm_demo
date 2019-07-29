% Load code and session variables
path(path, genpath('auto3dgm-matlab-gorgon'));
load('input/session.mat');

% Should have workspace variables ds, ga, pa, mst


% Generate globalized procrustes distances
k = 1;
proc_d     = zeros( ds.n , ds.n );
for ii = 1 : ds.n
    for jj = ii : ds.n
        if( ii == jj )
            continue;
        end
        [tmpR, proc_d( ii, jj)] = jprocrustes( ds.shape{ii}.X{k}*ga.P{ii} , ds.shape{jj}.X{k}*ga.P{jj} );
    end
end
proc_d = (proc_d+proc_d')/2;

% Generate 2D MST plot, using 2D MDS coordinates
plot_tree(proc_d, mst, ds.names, 'mds', ones(1,ds.n), 'MDS procrustes distances');
saveas(gcf, 'output/2D_MST_plot.png');

% Generate 3D MDS coordinates
coords = mdscale(proc_d,3)';

% Generate 3D MST plot, using 3D MDS coordinates
write_off_placed_shapes('output/map.off', coords, ds, ga, eye(3), mst);

% Tangent partial Procrustes coordinates PCA
tangent_pca(ds, ga, k);
saveas(gcf, 'output/tangent_procrustes_plot.png');