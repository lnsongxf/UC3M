'----------------------------------------------------------------------------------------
'Código escrito para el recuadro "Factores externos y fluctuaciones
'de la economía peruana"
'----------------------------------------------------------------------------------------
' (c) Carlos Rojas Quiroz
'----------------------------------------------------------------------------------------
' Lima, julio del 2019
'----------------------------------------------------------------------------------------
cd "C:\Users\User\Documents\21. Investigación\Choques externos\FKY Replica\Programa"
wfopen Base01.wf1
'----------------------------------------------------------------------------------------
'Características del sistema
'----------------------------------------------------------------------------------------
scalar nx			=4				'Número de variables externas
scalar nd    			=7				'Numero de variables domésticas
scalar nvar			=nx+nd   	'Numero de variables 
scalar nt          		=3           	'Términos determinísticos
scalar nparam 	=nx+nd+nt	'Número de parámetros estimados por ecuación
scalar numobs	=82 	     	'Número de observaciones
scalar lag			=1 			'Número de rezagos
scalar horz			=40        		'Longitud de la FIR
scalar nsimul   	=2000  		'Número de simulaciones del bootstrap
%initdate	="1998.1"
%enddate	="2018.2"
'----------------------------------------------------------------------------------------
'Variables del bloque exógeno
'----------------------------------------------------------------------------------------
series var1		=lgdpstar
series var2		=pic_star
series var3		=i_star 
series var4		=lipx
'----------------------------------------------------------------------------------------
'Variables del bloque doméstico
'----------------------------------------------------------------------------------------
series var5		=lgdp
series var6		=nonminning_sa
series var7		=minning_sa
series var8		=pic
series var9		=inom
series var10	=lrerbcrp
series var11	=ca_sa
'----------------------------------------------------------------------------------------
'Estimación del sistema
'----------------------------------------------------------------------------------------
smpl %initdate %enddate
include estimacion_sistema.prg
'Matriz de parámetros
coef F=C
matrix (nx+nd,nx+nd+nt) B1
scalar vv
for !m=1 to nx+nd
vv=(!m-1)*14
for !j=1 to nx+nd+nt
B1({!m},!j) = F(!j+vv)
next
next
'----------------------------------------------------------------------------------------
'Función Impulso Respuesta
'----------------------------------------------------------------------------------------
'Obtención de residuos
macro1.makeresids resid01 resid02 resid03 resid04 resid05 resid06 resid07 resid08 resid09 resid010 resid011
include construccion_fir.prg

for !k=1 to 11
matrix(horz,nvar) imp{!k}_f
imp{!k}_f=imp{!k}
next

'Boostrap para obtener Intervalo de Confianza
'include bootstrap.prg
