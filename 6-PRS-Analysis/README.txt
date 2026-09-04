## 6.1-PRS-Analysis-By-Age-Group.R

Input: phenotype data with PRS data
Output: results of analysis by age group
Function:
- fit linear model of metabolites ~ PRS using 1000 clustered bootstraps for each PRS in parallel
- extract relevant model metrics


## 6.2-PRS-Longtitudinal-Analysis.Returns

Input: phenotype data with PRS data
Output: results of Longtitudinal analysis
Function:
- fit linear mixed effects models, using natural cubic spline to model age, for each PRS in parallel
- extract relevant model metrics
