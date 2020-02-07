# loadshare
Prediction intervals for load sharing systems in accelerated life testing

This R package is an implementation of various prediction methods presented in the article 
    'Prediction intervals for load sharing systems in accelerated life testing' 
by K. Leckey, C.H. Müller, S. Szugat and R. Maurer (article submission FEB 2020)

The figures and tables in the article can be reproduced via the function calls

getFigure1(), getFigure2a(), getFigure2b(), getFigure3a(), getFigure3b(),
getFigure4a(), getFigure4b(), getSupFigure1a(), getSupFigure1b(), getSupFigure2a(),
getSupFigure2b(), getTable2()

Note that the run time of the functions varies depending on the background calculations 
that need to be made. Figures based on long-running simulations use saved .RData files
rather than repeating the simulations.
