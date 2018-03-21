# hazreg
This package implements novel models for survival analysis that all revolve around *hazard regularization*.
Four models are implemented. 
- The Age-Cohort Model is a well-known model in age-period-cohort analysis.
It is implemented with functions ending with `ac`.
- The Age-Cohort-Interaction Model is an extension of the age-cohort model.
It is implemented with functions ending with `aci`
- The Interaction Model is a generalization of the age-cohort model, in with no effect is estimated.
Only the hazard rate is estimated. It is implemented with functions ending with `interaction`.
- The Spatio-Temporal Piecewise Constant Hazard Model is a novel model aimed at inferring the hazard rate as 
a function of both time and a qualitative covariate representing the spatial location of individuals.

All these models have in common the same approach and the same regularization method at their code. 
This regularization is aimed at reducing overfitting (by reducing the variance) as well as improve the 
interpretability of the estimated hazard rate.
