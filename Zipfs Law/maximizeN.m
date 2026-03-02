function [n_max, best_array, best_slope, best_intercept, n_density] = maximizeN(ub, lb)
%Input metabolite concentration upper and lower bounds and determine
%maximized N after random permutations

% remove NaN
lb = lb(~isnan(lb),:);
ub = ub(~isnan(ub),:);

% Track the global best across all iterations
n_max = 0;
best_array = [];
iteration = NaN;
best_slope = 0;
best_intercept = 0;
n_density = zeros(100);

for k = 1:10000
% Generate random concentration list between 95% ci
rand_conc = (ub - lb).*rand(length(ub),1) + lb;
[rand_ordered, orig_indx] = sort(rand_conc, 'descend');
ranks_rand = find(rand_ordered);
rand_ordered = [ranks_rand, rand_ordered, orig_indx];


n = length(rand_ordered);

j = false;

slope     = NaN;
intercept = NaN;


while j == false && n>1
   top_N = rand_ordered(1:n,[1 2]);
   top_N_log = log10(top_N);
   mdl_log = fitlm(top_N_log(:,1),top_N_log(:,2));
   intercept = 10^(mdl_log.Coefficients{1,1}) ;
   slope = round(-mdl_log.Coefficients{2,1},2);

   % ----- dynamic bounds from the metabolite ranked #1 -----
   rank1_original_index = orig_indx(1);           % index in the *original* ub/lb
   rank1_lb = lb(rank1_original_index);
   rank1_ub = ub(rank1_original_index);
   % --------------------------------------------------------
 
   if (intercept <= rank1_ub) && (intercept >= rank1_lb) && (slope <= 1.1) && (slope >= 0.9)
       j = true;
       fprintf('Zipfs law fits the  data with N=%f\n', n);
       fprintf('The  slope is %f\n',slope);
       fprintf('The  intercept is %f\n', intercept);
       %slope_Jurkat = slope;
       %intercept_Jurkat = intercept;
       %N_Jurkat = n;
   else
       n = n - 1;
   end

end

n_density(k) = n;

if (n > n_max) || ...
   (n == n_max && abs(1 - slope) < abs(1 - best_slope))
    n_max      = n;
    best_array = rand_ordered; % this stores the original index
    iteration  = k;
    best_slope = slope;
    best_intercept = intercept;
else
end

end

end