function [avg_AbunPSNR, avg_AbunSSIM, avg_AbunSRE] = Evaluation(X_hat, XT, nl, nc, n, supp)

addpath("Utilities");

% 2D to 3D reshape
X_hat_im = reshape(X_hat',nl,nc,n); 
XT_im = reshape(XT',nl,nc,n);

% PSNR measure: average and individual abundance map
[avg_AbunPSNR, AbunPSNR_rec] = PSNR_cal(XT_im(:,:,supp), X_hat_im(:,:,supp), 0);

% SSIM measure: average and individual abundance map
AbunSSIM_rec = zeros(numel(supp),1);
for row = 1:numel(supp)
    AbunSSIM_rec(row,1) = SSIM(XT_im(:,:,supp(row)), X_hat_im(:,:,supp(row)));
end
avg_AbunSSIM = mean(AbunSSIM_rec);

% SRE measure: average and individual abundance map
AbunSRE_rec = zeros(numel(supp),1);
for row = 1:numel(supp)
    AbunSRE_rec(row,1) = 10*log10(norm(XT_im(:,:,supp(row)),'fro')^2 ...
                         / norm(X_hat_im(:,:,supp(row))-XT_im(:,:,supp(row)),'fro')^2);
end
avg_AbunSRE = mean(AbunSRE_rec);

end