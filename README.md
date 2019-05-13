# hazreg

This package implements models for survival analysis that revolve around *hazard regularization*.
Three models are implemented.

- The Age-Cohort-Interaction Model is an extension of the age-cohort model.
It is implemented by functions that end with `aci`

- The Interaction Model is a generalization of the age-cohort model, in which no effect is estimated.
Only the hazard rate is estimated. It is implemented by the functions that end with `interaction`.

- The Spatio-Temporal Piecewise Constant Hazard Model is a novel model aimed at inferring the hazard rate as 
a function of both time and a qualitative covariate representing the spatial location of individuals.

All these models have in common the same approach and the same regularization method at their core. 
This regularization is aimed at reducing overfitting as well as improve the 
interpretability of the estimated hazard rate.
