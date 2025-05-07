function outputs_analysis(mySect,faults,smears,faultDip,unit_plot,layerId,dim,Nsim,G0)
% 1 Visualize stratigraphy and fault (with thickness corresponding to 1st realization)
mySect.plotStrati(faults{1}.MatProps.thick, faultDip, unit_plot);  

% 2 Visualize intermediate variables
% We define a given parent material (id from 1 to n of materials in stratigraphy), 
% and generate histograms and correlation matrix plots.                                         
plotMatPropsHist(faults, smears, mySect, layerId, dim) 
% MatProps correlations
[R, P] = plotMatPropsCorr(faults, mySect, layerId, dim);

% 3 Visualize fault materials
% Visualization for one realization. Choice can be 'randm' (random), 'maxX' 
% (realization with maximum upscaled permeability in across the fault), 'minX', 
% 'maxZ' or 'minZ'.
% General fault materials and perm view
plotId = selectSimId('randm', faults, Nsim);                % simulation index
faults{plotId}.plotMaterials(mySect, G0) 

% 3.4. Visualize upscaled permeability
% Plot upscaled permeability distributions (all simulations)
plotUpscaledPerm(faults, dim)
end