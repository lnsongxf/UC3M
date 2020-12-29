'----------------------------------------------------------------------------------------
'Code for Metrics III project on the relationship between GDP and interest rates
'Country: Peru
'Authors: Tomasso D'Amelio and Carlos Rojas
'Date: January 2021
'----------------------------------------------------------------------------------------
cd "D:\Econometrics III"
import DataPeru.xlsx @freq m 2003.1
'----------------------------------------------------------------------------------------
'System's characteristics
'----------------------------------------------------------------------------------------
scalar k=2 'Number of endogenous variables
scalar dx=0 'Number of exogenous variables
scalar T=203 'Number of observations
scalar modelVar=0
%initdate="2003.1"
%enddate="2019.12"
'----------------------------------------------------------------------------------------
'Previous calculation
'----------------------------------------------------------------------------------------
series b=bc+bcX 'Total bank credit
series inom=iN*bc/b+ixN*bcX/b 'Average weighted interest rate
'Auxiliary variables
series var1=in-in(-1)
series var2=log(gdp/gdp(-1))*100
'----------------------------------------------------------------------------------------
'Optimal lag
'----------------------------------------------------------------------------------------
scalar maxlag=12
include optimallag
scalar maxlag=BICopt
'----------------------------------------------------------------------------------------
'VAR estimation
'----------------------------------------------------------------------------------------
smpl %initdate %enddate
!j=maxlag
if modelVar=0
equation eq01.ls var1 c var1(-1 to -!j) var2(-1 to -!j)
equation eq02.ls var2 c var1(-1 to -!j) var2(-1 to -!j)
else 
equation eq01.ls var1 c var1(-1 to -!j) var2(-1 to -!j) MA(1) MA(2) MA(3)
equation eq02.ls var2 c var1(-1 to -!j) var2(-1 to -!j) MA(1) MA(12)
next
'----------------------------------------------------------------------------------------
'Granger causality
'----------------------------------------------------------------------------------------
group Granger var1 var2
freeze(mode=overwrite, GrangerTest) Granger.cause(2)
delete Granger
'----------------------------------------------------------------------------------------
'Unit root
'----------------------------------------------------------------------------------------
freeze(mode=overwrite, DiffInterestRate_UR) var1.uroot(adf,none)
freeze(mode=overwrite, InterestRate_UR) in.uroot(adf,const)
freeze(mode=overwrite, DiffGDP_UR) var2.uroot(adf,const,trend)
freeze(mode=overwrite, GDP_UR) gdp.uroot(adf,const,trend)

