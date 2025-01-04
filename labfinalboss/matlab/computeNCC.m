function ncc = computeNCC(image1, image2)
    % Computes Normalized Cross-Correlation
    mean1 = mean(image1(:));
    mean2 = mean(image2(:));
    numerator = sum((image1(:) - mean1) .* (image2(:) - mean2));
    denominator = sqrt(sum((image1(:) - mean1).^2) * sum((image2(:) - mean2).^2));
    ncc = numerator / denominator;
end