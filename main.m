%%% 

[offres, asymmetry_phase] = OffResonance_Mapping(profiles, 4);

%%
[T1map, T2map, pdmap, b0map, amap, qmap, Mmap] = ...
  NLLS_Mapping(profiles, 5.7, 34, 30);

%%
selected_slice = 15;

subplot(1,2,1)
imagesc(squeeze(offres(:,:,selected_slice))), 
axis image, axis off, colormap hot, colorbar

subplot(1,2,2)
imagesc(squeeze(asymmetry_phase(:,:,selected_slice)),[-pi/6,pi/6]), 
axis image, axis off, colormap hot, colorbar

%%
figure
imagesc(squeeze(T1map(:,:,selected_slice)),[0 2000]), 
axis image, axis off, colormap inferno, colorbar
title("T1 map")

figure
imagesc(squeeze(T2map(:,:,selected_slice)),[0 150]), 
axis image, axis off, colormap inferno, colorbar
title("T2 map")

figure
imagesc(squeeze(pdmap(:,:,selected_slice)),[0 1e-1]), 
axis image, axis off, colormap inferno, colorbar
title("Proton Density Map")

figure
imagesc(squeeze(amap(:,:,selected_slice)),[0.8 0.99]), 
axis image, axis off, colormap copper, colorbar
title("Elliptical parameter a")

figure
imagesc(squeeze(qmap(:,:,selected_slice)),[0.4 0.8]), 
axis image, axis off, colormap copper, colorbar
title("Elliptical parameter q")

figure
imagesc(squeeze(Mmap(:,:,selected_slice)),[0 1e-2]), 
axis image, axis off, colormap copper, colorbar
title("Elliptical parameter M")




