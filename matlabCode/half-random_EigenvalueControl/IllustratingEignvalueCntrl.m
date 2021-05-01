% Little script illustrating eigenvalue control (related to Fig. 1e, Fig. S2,
% and the 'half-random' results for the correlations between eigenvectors
% in Fig. S2 and Fig. 2).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Creating the cortico-cortical weight matrix %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setting the seed
s=1;
rng(s);

% Size of the network
N = 500;

% Original matrix, and its right eigenvectors and eigenvalues
M = randn(N, N)/sqrt(N);    
[R, D] = eig(M); % R: original right eigenvectors
E = diag(D);   % E: original eigenvalues
L=inv(R);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Reproducing Fig. 1e %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of controlled eigenvalues
n = 20;
% Position of the new desired eigenvalues
Ei = 1 + 1*1i*linspace(-1, 1, n);

%%%%%%% Creating a half-random rank-one perturbation u*v.' to impose the
%%%%%%% desired eigenvalues

% Creating the matrix P (from Eqs. 2 and 10)
P = zeros(n, N);
for j=1:n
    P(j, :) = 1./(Ei(j) - E);
end

% Create the vector d from Eq. 12, which fixes all of the eigenvalues
d = pinv(P)*ones(n, 1);

% Choosing a random u
u=randn(N,1);

% Computing the corresponding vector v by using Eq. 2 (see also Eq. 12)
a=L*u;
b=(d./a);
ad_vT=(b.')*L;


%%%%%%% Computing the spectrum of the new matrix
Mp = M + u*ad_vT;
[Rp, Dp] = eig(Mp);
Ep = diag(Dp);

%%%%%%% Ploting
figure(1)
plot(real(E), imag(E), 'r.','MarkerSize',8);   % plot old eigenvalues
hold on;
plot(real(Ep), imag(Ep), 'o','Color', [0.5,0,0.5],'MarkerSize',6); % plot new eigenvalues
plot(real(Ei), imag(Ei), 'g+','MarkerSize',3); % plot desired eigenvalues
hold off;
xlim([-1.2,1.2])
ylim([-1.2,1.2])
axis equal
xlabel('Real part of eigenvalues')
ylabel('Imaginary part of eigenvalues')
title('Eigenspectrum')
legend({'Initial cortical','Perturbed','Target'},'Location','northwest')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Showing that it is equally possible to control u or v %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% As long as d is fixed, this fixes the whole eigenspectrum, and
% it is always possible to fix d by using an arbitrary u  - or an
% arbitrary v - and adjusting the other vector of the loop
% correspondingly.

% We keep the same d as before, but choose a random v and will
% adapt u accordingly.

new_vT=randn(1,N);

% Computing the corresponding vector u (see Eq. 12: a=L*u ; b=(v.')*R )
b=new_vT*R;
a = d./(b.'); 

ad_u = R*a;

%%%%%%% Computing the spectrum of the new matrix
Mp2 = M + ad_u*new_vT;
[Rp2, Dp2] = eig(Mp2);
Ep2 = diag(Dp2);

% Show that the matrices Mp and Mp2 are different
%max(abs(Mp(:)-Mp2(:)))./mean(abs(Mp(:)))
% 16.4308 -> the matrices are different

%%%%%%% Ploting
figure(2)
plot(real(E), imag(E), 'r.','MarkerSize',8);   % plot old eigenvalues
hold on;
plot(real(Ep), imag(Ep), 'o','Color', [0.5,0,0.5],'MarkerSize',6); % plot new eigenvalues of Mp
plot(real(Ep2), imag(Ep2), 'o','Color', [0,0,0],'MarkerSize',3); % plot new eigenvalues of Mp2
plot(real(Ei), imag(Ei), 'g+','MarkerSize',3); % plot desired eigenvalues
hold off;
xlim([-1.2,1.2])
ylim([-1.2,1.2])
axis equal
xlabel('Real part of eigenvalues')
ylabel('Imaginary part of eigenvalues')
title('Eigenspectrum')
legend({'Initial cortical','Perturbed 1','Perturbed 2','Target'},'Location','northwest')
