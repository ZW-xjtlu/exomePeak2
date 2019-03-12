library(exomePeak2)
library(Rsamtools)

n = 100
Rn = 0.645
A = round( 1+(qnorm(.025)^2)/n, 4)
B = round(-(2*Rn + (qnorm(.025)^2)/n), 4)
C = round(Rn^2, 4)

(-B + sqrt(B^2 - 4*A*C))/(2*A)

(-B - sqrt(B^2 - 4*A*C))/(2*A)

