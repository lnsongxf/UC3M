%initdate="2003.1"
%enddate="2019.12"
%initdate2 = @otod(@dtoo(%initdate)+maxlag)
smpl %initdate2 %enddate
for !j=1 to maxlag
scalar P=!j
'VAR model
if modelVar==0 then
equation eq01.ls var1 c var1(-1 to -!j) var2(-1 to -!j)
equation eq02.ls var2 c var1(-1 to -!j) var2(-1 to -!j)
else
equation eq01.ls var1 c var1(-1 to -!j) var2(-1 to -!j) MA(1) MA(2) MA(3)
equation eq02.ls var2 c var1(-1 to -!j) var2(-1 to -!j) MA(1) MA(12)
endif
' Working on residuals
eq01.makeresids resid01 
eq02.makeresids resid02
stom(resid01,e01)
stom(resid02,e02)
matrix(2,@rows(e01)) ee
matplace(ee,@transpose(e01),1,1)
matplace(ee,@transpose(e02),2,1)
matrix(2,2) sume
sume=ee*@transpose(ee)
delete resid01 resid02 e01 e02 ee
' Residual covariance: ee*ee'
scalar sume2=@det(1/T*sume)
' Log-likelihood
scalar loglike=-T/2*(k*(1+log(2*@pi))+log(sume2))
' Information criteria
matrix(maxlag,1) akaike
matrix(maxlag,1) schwartz
matrix(maxlag,1) hannanquinn
akaike(!j,1)=-2*loglike/T+2*(k*(p*k+dx))/T
schwartz(!j,1)=-2*loglike/T+(k*(p*k+dx))*log(T)/T
hannanquinn(!j,1)=-2*loglike/T+2*(k*(p*k+dx))*log(log(T))/T
next
scalar AICOpt = @cimin(akaike)
scalar BICOpt = @cimin(schwartz)
scalar HQOpt = @cimin(hannanquinn)


