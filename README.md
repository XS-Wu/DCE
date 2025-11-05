# Copyright © 2025 Xinsheng Wu, Fudan University


# Overview
This repository contains R code to generate D-optimal discrete choice experiment (DCE) designs for PrEP preference studies using the idefix package.
It builds 60 tasks with 2 alternatives (60×2), decodes human-readable questionnaires, optionally splits tasks into balanced blocks, and exports diagnostics and helper tables.

# Requirements
- R (≥ 4.2 recommended)
- Packages: idefix (script calls library(idefix))

# Quick Start
1) Install dependency
   install.packages("idefix")
2) Run the script (RStudio or command line)
   source("DCEscript.R") 
3) Main outputs (matching the script naming):
   - DCE_60x2_idefix_questionnaire_A.csv
   - DCE_60x2_idefix_designMatrix_A.csv
   - Blocks_A_Block01.csv, Blocks_A_Block02.csv, Blocks_A_Block03.csv
   - Same for label B
   - Table_Attributes_and_Levels.csv (attribute–level table)
   - Table_Priors.csv (prior vector)


# License
All rights reserved. See the copyright notice above. For usage beyond personal research, please contact the author for permission.

