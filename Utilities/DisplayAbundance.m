function DisplayAbundance(X_hat, XT, nl, nc, n, p, supp, Method)

X_hat_im = reshape(X_hat',nl,nc,n); 
XT_im = reshape(XT',nl,nc,n);
        
fig = figure();

[ha, pos] = tight_subplot(p,2,[.025 .015],[.05 .05],[.05 .15]);

for j = 1:p
    axes(ha(1+(j-1)*2));
    imagesc(XT_im(:,:,supp(j)), [0 1])
    axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])

    axes(ha(2+(j-1)*2));
    imagesc(X_hat_im(:,:,supp(j)), [0 1])
    axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])
  
    originalSize2 = get(gca, 'Position');
    h = colorbar; 
    set(ha(2+(j-1)*2), 'Position', originalSize2);
    set(h,'fontsize',8);    
end
colormap jet

numAlgs = 2;
for j = 1:p
    axes(ha((j-1)*numAlgs + 1)); ylabel(sprintf('Abundance %i',j),'interpreter','latex');
end

axes(ha(end-1)); xlabel('GT','Interpreter','latex')
axes(ha(end)); xlabel(Method,'Interpreter','latex')

end