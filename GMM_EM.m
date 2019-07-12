% GMM EM

function [] = GMM_EM()

    K = 30;

    Q = zeros(21604, K);
    t = rand(K, 1);
    w = t ./ sum(t);
    mu = rand(K, 2) .* 100;
    E = zeros(21604, K);
    X_mat = repmat(c1__, K);
    
    for i = 1:100
    
        % E-step
        mu_mat = repmat();
        
        % M-step
        
    end
end
