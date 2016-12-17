#!/bin/bash

echo "Running exponential_cohort_model.R..."
rscript exponential_cohort_model.R
echo "Running linear_cohort_model.R..."
rscript linear_cohort_model.R
echo "Running multi_grade_model.R..."
rscript multi_grade_model.R
echo "Running multi_grade_model_interactions.R..."
rscript multi_grade_model_interactions.R
echo "Running multi_grade_model_quadratic.R..."
rscript multi_grade_model_quadratic.R
