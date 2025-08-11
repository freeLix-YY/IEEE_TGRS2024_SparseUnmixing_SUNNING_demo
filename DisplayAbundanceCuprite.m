function DisplayAbundanceCuprite(X_hat, nl, nc, n, Method)

material_idx = [420, 336, 297];
range_mat = [0 0.38; ... % Alunite
             0 0.30; ... % Buddingtonite
             0 0.65];    % Chalcedony

X_hat_im = reshape(X_hat',nl,nc,n); 
        
fig = figure();

[ha, pos] = tight_subplot(3,1,[.025 .015],[.05 .05],[.05 .15]);

for j = 1:numel(material_idx)
    axes(ha(j));
    imagesc(X_hat_im(:,:,material_idx(j)), range_mat(j,:))
    axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])
  
    originalSize2 = get(gca, 'Position');
    h = colorbar; 
    set(ha(j), 'Position', originalSize2);
    set(h,'fontsize',8);    
end
colormap jet

for j = 1:numel(material_idx)
    axes(ha(j)); ylabel(sprintf('Abundance %i',j),'interpreter','latex');
end

axes(ha(end)); xlabel(Method,'Interpreter','latex')

end