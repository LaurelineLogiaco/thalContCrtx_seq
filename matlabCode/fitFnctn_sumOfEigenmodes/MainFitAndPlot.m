% Script to find the parameters for the combination of a small number of
% basis functions into a desired network output.

vec_Halfn=[2,3,4,5,10];
% vectors with several numbers of eigenmodes (several fits with increasing precision)

cmap = colormap(parula(length(vec_Halfn)));
% different colors for the different fits with increasing precision

n_diff_randInitSeeds=100;
% We will start each fit from 100 random initial conditions

nonlincs = @(x) constraints(x,3,0.05,2);
% The constraints when fitting the desired target function:
% 1) in order to have a robust fit, the complex amplitudes of the basis
% functions cannot be too large (which can become an issue when using a
% relatively large number of basis functions). Here, 3 is the maximum norm
% of the vector of complex amplitudes.
% 2) for the perturbed network to be able to be extended to include the
% desired eigenvalues corresponding to the timescales of the basis
% functions that can fit the desired output well, these eigenvalues can't
% be too close to one another (which would create strong eigenvector
% correlations - here we should the conservative minimum distance between
% eigenvalues of 0.05) or too far away from the center of an arbitrary
% initial eigenspectrum (which is typically a spheroid, so here we just
% impose a maximal distance between the desired eigenvalues of 2). 

dt = 0.1;
% Timestep.

A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
% The 'default' constraints for the optimization of the are empty because
% all of our constraints are customized through the function nonlincs.




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Simple oscillatory pattern %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load oscill_patterns

f=patterns{4,1};
% One example pattern

T=length(f);

t = dt * (1:T);
% vector of times

vect_loss_fun_SeanPatt=zeros(length(vec_Halfn),1);
% vector with the final loss for all numbers of basis functions used

cell_opt_sols_SeanPatt=cell(length(vec_Halfn),1);
% cell with the optimized parameters for all numbers of basis functions

mat_fvs_SeanPatt=zeros(length(vec_Halfn),T);
% The weighted sum of basis function that approximates f for different
% numbers of basis functions used.

for l=1:length(vec_Halfn)
    
    s=1;
    rng(s);
    % fixing the seed
    
    n=vec_Halfn(l);
    
    ub = [[zeros(n,1);inf*ones(n,1)], inf*ones(2*n,1)];
    % Upper bounds:
    % zeros for the timescales of the basis functions (we still got good results while ensuring more stable dynamics).
    % infinity for all other parameters (though the amplitudes will be constrained with nonlincs)
    
    options = optimoptions('fmincon','MaxFunctionEvaluations',1000*10*n,'MaxIterations',1000*10*n,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',false,'Display','final-detailed','CheckGradients',false,'FiniteDifferenceType','central');
    
    vec_interm_valTrgt=zeros(1,n_diff_randInitSeeds);
    % reminding all losses for all different random starts of the optimization
    
    cell_interm_xOpts=cell(1,n_diff_randInitSeeds);
    % reminding all parameters for all different random starts of the optimization
    
    % loops over the different initial random seeds
    for k=1:n_diff_randInitSeeds
        
        % while statement to make sure that the parameters only include
        % negative exponential rates.
        
        while true
            
            % random initialization
            x0 = randn(2 * n, 2);
            
            % nonlinear constrained optimization
            x = fmincon(@(x) cost_func(x,f,dt),x0,A,b,Aeq,beq,lb,ub,nonlincs,options);
            
            % verifying that there are no positive exponential rates, if
            % not, keep optimizing.
            tmp=x(1:n,1)+1i*x(n+1:end,1);
            if any(real(tmp)>=0)
                continue
            end
            break;
        end

        vec_interm_valTrgt(k)=cost_func(x,f,dt);
        % remind the final loss after optimization with these initial
        % values of the parameters.
        
        cell_interm_xOpts{k}=x;
        % remind the final parameters that the algorithm has converged to
        % with these initial parameter values.
    end
    
    [min_val,min_ind]=min(vec_interm_valTrgt);
    % find the best optimization among all initial parameter values.

    f_tar = EignMd_Approx(cell_interm_xOpts{min_ind},T,dt);
    % compute the best function for the best parameters found given the
    % current number of basis functions that are combined, and store it
    % below.
    mat_fvs_SeanPatt(l,:)=f_tar;
    
    vect_loss_fun_SeanPatt(l)=sum((f_tar-f).^2);
    % remind the smallest loss found given the current number of basis
    % functions that are combined.
    
    cell_opt_sols_SeanPatt{l}=cell_interm_xOpts{min_ind};
    % remind the parameters for the smallest loss found.
    
end

% uncomment the three lines below to save the parameters

%save('cell_opt_sols_SeanPatt.mat','cell_opt_sols_SeanPatt')
%save('vect_loss_fun_SeanPatt.mat','vect_loss_fun_SeanPatt')
%save('mat_fvs_SeanPatt.mat','mat_fvs_SeanPatt')


%%%%%%%%%%%%%%%%%%%%%%
% plotting the results
%%%%%%%%%%%%%%%%%%%%%%

% plotting the target function vs. the fit

figure()
hold on
plot(t,f,'linewidth',3,'color','k')
% plotting the target function

% loop over the number of basis functions
for i=length(vec_Halfn):-1:1
    plot(t,mat_fvs_SeanPatt(i,:),'linewidth',1,'color',cmap(i,:))
    % plotting the best found fit for different numbers of basis functions
end
%saveas(gcf,'temporalFunctions_SeanPatt.fig')
xlabel('Time')
title({'Compa. smooth target (black)'; 'and fits with diff. # of basis functions'})

% plotting the best eigenvalues found (related to the timescales of the
% basis functions)
figure()
hold on
for i=length(vec_Halfn):-1:1
    n=vec_Halfn(i);
    plot([cell_opt_sols_SeanPatt{i}(1:n,1);cell_opt_sols_SeanPatt{i}(1:n,1)],[cell_opt_sols_SeanPatt{i}((n+1):end,1);-cell_opt_sols_SeanPatt{i}((n+1):end,1)],'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',cmap(i,:))
end
%saveas(gcf,'eignvals_SeanPatt.fig')
xlabel('Real part of eigenvalues')
ylabel('Imaginary part of eigenvalues')
title({'Eigenvalues of the basis functions that fit the' ; 'smooth target (for diff. # of basis functions)'})


% plotting the complex weights for the basis functions that leads to the
% best found fit of the desired target.
figure()
hold on
for i=length(vec_Halfn):-1:1
    n=vec_Halfn(i);
    plot([cell_opt_sols_SeanPatt{i}(1:n,2);cell_opt_sols_SeanPatt{i}(1:n,2)],[cell_opt_sols_SeanPatt{i}((n+1):end,2);-cell_opt_sols_SeanPatt{i}((n+1):end,2)],'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',cmap(i,:))
end
%saveas(gcf,'Amp_SeanPatt.fig')
xlabel('Real part of the weights')
ylabel('Imaginary part of the weights')
title({'Weights of the basis functions that fit the' ; 'smooth target (for diff. # of basis functions)'})



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Hard stepwise pattern %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load cell_stpd_trajectX.mat


f=cell_stpd_trajectX{10};

T=length(f);
t = dt * (1:T);

vect_loss_fun_SqrPatt1=zeros(length(vec_Halfn),1);

cell_opt_sols_fun_SqrPatt1=cell(length(vec_Halfn),1);

mat_fvs_SqrPatt1=zeros(length(vec_Halfn),T);

for l=1:length(vec_Halfn)
    
    s=1;
    rng(s);
    
    n=vec_Halfn(l);
    
    ub = [[zeros(n,1);inf*ones(n,1)], inf*ones(2*n,1)];
    options = optimoptions('fmincon','MaxFunctionEvaluations',1000*10*n,'MaxIterations',1000*10*n,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',false,'Display','final-detailed','CheckGradients',false,'FiniteDifferenceType','central');
    
    vec_interm_valTrgt=zeros(1,n_diff_randInitSeeds);
    cell_interm_xOpts=cell(1,n_diff_randInitSeeds);
    
    for k=1:n_diff_randInitSeeds

        while true
            x0 = randn(2 * n, 2);
            x = fmincon(@(x) cost_func(x,f,dt),x0,A,b,Aeq,beq,lb,ub,nonlincs,options);
            tmp=x(1:n,1)+1i*x(n+1:end,1);
            if any(real(tmp)>=0)
                continue
            end
            break;
        end

        vec_interm_valTrgt(k)=cost_func(x,f,dt);
        cell_interm_xOpts{k}=x;
    end
    
    [min_val,min_ind]=min(vec_interm_valTrgt);

    f_tar = EignMd_Approx(cell_interm_xOpts{min_ind},T,dt);
    
    vect_loss_fun_SqrPatt1(l)=sum((f_tar-f).^2);
    
    cell_opt_sols_fun_SqrPatt1{l}=cell_interm_xOpts{min_ind};

    mat_fvs_SqrPatt1(l,:)=f_tar;
end

%save('cell_opt_sols_fun_SqrPatt1.mat','cell_opt_sols_fun_SqrPatt1')
%save('vect_loss_fun_SqrPatt1.mat','vect_loss_fun_SqrPatt1')
%save('mat_fvs_SqrPatt1.mat','mat_fvs_SqrPatt1')

%%%%%%%%%%%%%%%%%%%%%%
% plotting the results
%%%%%%%%%%%%%%%%%%%%%%


figure()
hold on
plot(t,f,'linewidth',3,'color','k')
for i=length(vec_Halfn):-1:1
    plot(t,mat_fvs_SqrPatt1(i,:),'linewidth',1,'color',cmap(i,:))
end
%saveas(gcf,'temporalFunctions_SqrPatt1.fig')
xlabel('Time')
title({'Compa. step target (black)'; 'and fits with diff. # of basis functions'})

figure()
hold on
for i=length(vec_Halfn):-1:1
    n=vec_Halfn(i);
    plot([cell_opt_sols_fun_SqrPatt1{i}(1:n,1);cell_opt_sols_fun_SqrPatt1{i}(1:n,1)],[cell_opt_sols_fun_SqrPatt1{i}((n+1):end,1);-cell_opt_sols_fun_SqrPatt1{i}((n+1):end,1)],'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',cmap(i,:))
end
%saveas(gcf,'eignvals_SqrPatt1.fig')
xlabel('Real part of eigenvalues')
ylabel('Imaginary part of eigenvalues')
title({'Eigenvalues of the basis functions that fit the' ; 'step target (for diff. # of basis functions)'})

figure()
hold on
for i=length(vec_Halfn):-1:1
    n=vec_Halfn(i);
    plot([cell_opt_sols_fun_SqrPatt1{i}(1:n,2);cell_opt_sols_fun_SqrPatt1{i}(1:n,2)],[cell_opt_sols_fun_SqrPatt1{i}((n+1):end,2);-cell_opt_sols_fun_SqrPatt1{i}((n+1):end,2)],'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',cmap(i,:))
end
%saveas(gcf,'Amp_SqrPatt1.fig')
xlabel('Real part of the weights')
ylabel('Imaginary part of the weights')
title({'Weights of the basis functions that fit the' ; 'step target (for diff. # of basis functions)'})
