function ngc = computeNGC(image1, image2)
    % Computes Normalized Gradient Correlation
    [gx1, gy1, gz1] = gradient(double(image1));
    [gx2, gy2, gz2] = gradient(double(image2));
    dotProduct = sum((gx1(:) .* gx2(:)) + (gy1(:) .* gy2(:)) + (gz1(:) .* gz2(:)));
    norm1 = sqrt(sum(gx1(:).^2 + gy1(:).^2 + gz1(:).^2));
    norm2 = sqrt(sum(gx2(:).^2 + gy2(:).^2 + gz2(:).^2));
    ngc = dotProduct / (norm1 * norm2);
end
