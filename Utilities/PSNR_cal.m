%% PSNR calculation
function [cpsnr,psnr] = PSNR_cal(org,recon,skip)
    %skip : to skip some boundary lines
    org=org(skip+1:end-skip,skip+1:end-skip,:);
    recon=recon(skip+1:end-skip,skip+1:end-skip,:);
      [m, n,~]=size(org);
    
    if strcmp(class(org),class(recon))
        sse=squeeze(sum(sum((org-recon).^2))); %square sum of error  
        mse=sse./(m*n);  %mean square error of each band.
        maxval=squeeze(max(max(org)));
        psnr= 10*log10( (maxval.^2) ./mse);
        cpsnr=mean(psnr);
    end
end