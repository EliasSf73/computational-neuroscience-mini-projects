
% Prob1d: numerical vs analytic entropy for exponential (geometric) spike counts

  % 1) parameters
  r      = 15;         
  T      = 1;          
  N      = 100;        
  lambda = r * T;      

  % 2) sample from discrete exponential (geometric) distribution
  rng(45);
  p_geo  = 1/(1+lambda);         
  counts = geornd(p_geo, [N,1]);  

  % 3) build empirical pmf
  edges      = -0.5 : 1 : (max(counts)+0.5);
  count_hist = histcounts(counts, edges);
  p_emp      = count_hist / N;  
  n_vals     = 0:(numel(p_emp)-1);

  % 4) numerical entropy
  idx   = p_emp > 0;
  S_num = -sum( p_emp(idx) .* log2(p_emp(idx)) );

  % 5) analytic entropy
  S_an  = log2(1+lambda) + lambda*log2(1 + 1/lambda);

  % 6) absolute error
  delta = abs(S_num - S_an);

  % display
  fprintf('Exponential code (lambda=%.1f, N=%d)\n', lambda, N);
  fprintf('  S_num = %.3f  bits\n', S_num);
  fprintf('  S_an  = %.3f  bits\n', S_an);
  fprintf('  |ΔS|  = %.3f  bits\n', delta);

