function expect_rmse=uOpt_ToMin_NoisyDev(u,cell_params)

% Computes the cost in Eq. 22, keeping the eigenvalues of the
% motif-specific dynamics fixed (by fixing d from Eq. 13) but changing the
% eigenvectors to stabilize the readout direction. Note that, as we
% explained below equation 9, the choice of u vs. v as remaining degrees of
% freedom after fixing d, that can then be optimized in this function to
% improve readout robustness, is arbitrary. Here, we optimize over u.

InitRightEigVct=cell_params{1};
% matrix of right eigenvectors of the unperturbed matrix

readout=cell_params{2};
% readout vector

AmpPhase0=cell_params{3};
% complex amplitudes for the basis functions

t_end=cell_params{4};
% duration of the motif

pct_max=cell_params{5};
% just a scalar such that the output of the cost function is expressed as
% the expected deviation when the initial noise scales as a certain
% percentage of the magnitude of the motif itself.

dim=cell_params{6};
% the size of the cortical network (the nber of rows and columns of
% InitRightEigVct, and the number of elements of readout)

mat_timescale=cell_params{7};
% corresponds to the matrix \Lambda in Eq. 22

mat_Q_dfNewOldEgnvl_smDfIds=cell_params{8};
% to quickly compute the eigenvectors of the perturbed matrix from the
% eigenvectors of the isolated cortical matrix, without needing an
% expensive eigendecomposition of the perturbed matrix.
% Specifically, we compute the un-normalized version of Eq. 18, at the top of the
% 22nd page of the manuscript (vectors indexed with an
% open circle on top in the manuscript). This involves, for the ith right eigenvector, the
% ith column of the matrix A multiplied by the ith entry of the diagonal of
% the matrix A (corresponding to the terms including the eigenvalues in the un-normalizd version of the eigenvectors).
% This gives a column vector for each right eigenvectors. The matrix Q is
% the concatenation of all these right eigenvectors.

InitLeftEigVct=cell_params{9};
% matrix of left eigenvectors of the unperturbed matrix

LeftInitEigVc_Mult_u=InitLeftEigVct*u;
% projection of u on the left initial eigenvactors

vct_weights_NewRightEig_intoOld=(repmat(LeftInitEigVc_Mult_u,1,dim)./repmat(LeftInitEigVc_Mult_u.',dim,1)).*mat_Q_dfNewOldEgnvl_smDfIds;
% these weights correspond to how to form the right eigenvector of the full
% (perturbed) matrix from the right eigenvectors of the original matrix -
% this is the unnormalized version of the eigenvectors, indexed with an
% open circle on top in the manuscript.

FinRightEigVct=InitRightEigVct*vct_weights_NewRightEig_intoOld;
% This is the unnormalized matrix of right eigenvectors for the perturbed
% matrix.

% For stability of the inversion, we will normalize the eigenvectors.
FinRightEigVct=FinRightEigVct./repmat(sqrt(sum(FinRightEigVct.*conj(FinRightEigVct))),dim,1);

% weights to recover the readout as a linear combination of the left
% eigenvectors of the perturbed matrix.
wgts_LEgnvct_inRdt=readout*FinRightEigVct;

% initial complex amplitude (equation 16)
A0=FinRightEigVct*diag(1./wgts_LEgnvct_inRdt)*AmpPhase0;

% Left eigenvectors (i.e., inv(FinRightEigVct)) times A0
vct_FinLeftEigVctTiA0=FinRightEigVct\A0;

% Matrix with scalar product between right eigenvectors, entering in the
% equation of the noise (Eq. 23)
fRt_fR=(FinRightEigVct.')*FinRightEigVct;

% This is equation 23, except un-normalized (no division by the number of
% units or the duration of the activity transient).
IntNrmActs=(vct_FinLeftEigVctTiA0.')*(fRt_fR.*mat_timescale)*vct_FinLeftEigVctTiA0;

%%%% Finally: Amplitude of noise - equivalent of the square root of Eq. 23,
%%%% except scaled by a factor sqrt(dim).
AmpNoise=pct_max*sqrt((IntNrmActs/t_end));
% pct_max is superfluous but it allows the output numbers to be in a
% meaningful scale, and it is a negligible 'weight'.

% The matrix of scalar product between the left eigenvectors (can be
% obtained by directly inverting the matrix of scalar product between right
% eigenvectors, as inv(A B) = inv(B) inv(A).
fL_fLt=inv(fRt_fR);

% Finally, the square root of Eq. 22, rescaled
expect_rmse = real((AmpNoise/sqrt(dim)) *sqrt(wgts_LEgnvct_inRdt * ((fL_fLt.*mat_timescale) * wgts_LEgnvct_inRdt.')/t_end));
