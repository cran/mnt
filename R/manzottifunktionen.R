spere=function(x){
  10*(x[,1]^3*x[,2]^3+x[,1]^3*x[,3]^3+x[,2]^3*x[,3]^3)-3*(x[,1]*x[,2]^5+x[,1]*x[,3]^5+x[,2]*x[,3]^5+x[,2]*x[,1]^5+x[,3]*x[,1]^5+x[,3]*x[,2]^5)
}

sphere2=function(x){
  d=length(x)
  A=5*x^3%*%t(as.matrix(x^3))
  diag(A)=rep(0,d)
  B=-3*x^5%*%t(as.matrix(x))
  diag(B)=rep(0,d)
  return(sum(A)+sum(B))
}


sphere3=function(data){
  n=dim(data)[1]
  ret=rep(0,n)
  for (j in 1:n) ret[j]=sphere2(data[j,])
  return(ret)
}



##Normmierungskonstante
##input: d (Dimension) q (Grad)
##output: nat?rliche Zahl d(d+2)(d+4)...
norm=function(d,q){
  x=c(0:(q-1))*2+d  #Vector der L?nge q mit Zahlen d,d+2,d+4,...
  ret=prod(x)		#Produkt
  return(ret)
}


###############################################################################################################
## Berechnet die skalierten Residuen
## input: n x d Daten-Matrix n(Stichprobenumfang) d (Dimension)
## output: n x d Matrix
normdaten=function(data){
  n=dim(data)[1] 						# Stichprobenumfang
  meanary=t(array(colMeans(data),dim(t(data))))   # Mittelwertmatrix
  zdata=data-meanary					# zentrierte Stichprobe
  Sn=stats::cov(data)*(n-1)/n					# empirische Kovarianzmatrix
  invSn=solve(Sn)						# inverse empirische Kovarianzmatrix
  V=eigen(invSn)$vectors					# Eigenwerte Matrix
  e=eigen(invSn)$values					# Eigenwerte
  sqrtSn=V%*%diag(sqrt(e))%*%t(V)			# positive quadratwurzel aus invSn
  ret=sqrtSn %*% t(zdata)					# normierte Stichprobe
  return(t(ret))
}

###################################################################################################################
## Berechnung der Winkelkomponenten U
## input: n x d Matrix
## output: n x d Matrix

U=function(data){
  Prod=data*data						# Quadrat komponentenweise
  norm=array(sqrt(rowSums(Prod)),dim(data))		# Norm-Matrix
  ret=data/norm						# normierte Daten
  return(ret)
}

#######################################################################################################################
## Berechnung der radialen komponente R
## input: n x d Matrix
## output: n x d Matrix

R=function(data){
  Prod=data*data			# Quadrat komponentenweise
  ret=sqrt(rowSums(Prod))		# Norm-Vektor
  return(ret)
}

## harmonische Kugelfunktionen von Grad 2 bis 4

f2=function(x,d,A){ 	X=t(t(x))%*%x
				n=sqrt(norm(d,2))
				R=X[A]*n
				return(R)
}

f3=function(x,d,H){ 	y=t(x^2)
				Y=H%*%t(y)
				r=c(2:d)
				r=sqrt(norm(d,2)/(2*r*(r-1)))
				R=Y[2:d]*r
				return(R)
}

f4=function(x,d,A){ 	E=array(0,dim=c(d,d,d))
				X=t(t(x))%*%x
				x2=x*sqrt(norm(d,3))
				for (i in 1:d){
					E[ , ,i]=X*x2[i]
				}
				R=E[A]
				return(R)
}



f5=function(x,d,A,H){ 	y=t(x^2)
				y=H%*%t(y)
				r=c(1:d)
				r=sqrt(norm(d,3)/(2*r*(r-1)))
				y=y*r
				X=y%*%x
				R=X[A]
				return(R)
}


f6=function(x,d,A,H){ 	y=t(x^2)
				y=H%*%t(y)
				r=c(1:d)
				r=sqrt(norm(d,3)/(2*(r+2)*(r+1)))
				y=y*r
				X=t(t(x))%*%t(y)
				R=X[A]
				return(R)
}


f7=function(x,d,H){ 	y3=t(x^3)
				y2=t(x^2)
				y2=H%*%t(y2)
				r=c(1:d)
				Y=(r-1)*y3-3*t(x)*t(y2)
				r2=sqrt(norm(d,3)/(6*(r+2)*(r-1)))
				R=Y[2:d]*r2[2:d]
				return(R)
}

f8=function(x,d,H){ 	y4=t(t(x^4))
				y2=t(t(x^2))
				s=H%*%y2
				s2=s^2
				r=c(1:d)
				Y=6*(r+1)*y2*s-(r^2-1)*y4-3*s2
				r2=sqrt(norm(d,4)/(24*(r-1)*(r+1)*(r+2)*(r+4)))
				R=Y[2:d]*r2[2:d]
				return(R)
}

f9=function(x,d,A,H){ 	y=t(x^2)
				y=H%*%t(y)
				y=y*x
				r=c(1:d)
				r=sqrt(norm(d,4)/(6*(r+1)*(r+4)))
				y=y*r
				X=t(t(x))%*%t(y)
				R=X[A]
				return(R)
}

f10=function(x,d,A,H){ 	E=array(0,dim=c(d,d,d))
				X=t(t(x))%*%x
				y2=t(t(x^2))
				y=H%*%y2
				r=c(1:d)
				r=sqrt(norm(d,4)/(2*(r+3)*(r+4)))
				y=y*r
				for (i in 1:d){
					E[ , ,i]=X*y[i]
				}
				R=E[A]
				return(R)
}

f11=function(x,d,A,H1,H2){ 	y2=t(x^2)
					s=H1%*%t(y2)
					y=H2%*%t(y2)
					r=c(1:d)
					r2=sqrt(1/(4*r*(r-1)))
					r3=sqrt(1/((r+3)*(r+4)))
					s=s*r3
					y=y*r2
					X=y%*%t(s)
					X=X*sqrt(norm(d,4))
					R=X[A]
					return(R)
}

f12=function(x,d,A){ 	E=array(0,dim=c(d,d,d,d))
				X=t(t(x))%*%x
				x2=x*sqrt(norm(d,4))
				for (i in 1:d){
					for(j in 1:d){
						E[ , ,i,j]=X*x[i]*x2[j]
					}
				}
				R=E[A]
				return(R)
}


f13=function(x,d,A,H){ 	E=array(0,dim=c(d,d,d))
				y2=t(t(x^2))
				y=H%*%y2
				r=c(1:d)
				r=sqrt(norm(d,4)/(2*r*(r-1)))
				y=y*r
				X=y%*%x
				for (i in 1:d){
					E[ , ,i]=X*x[i]
				}
				R=E[A]
				return(R)
}

f14=function(x,d,A,H){ 	E=array(0,dim=c(d,d,d))
				y2=t(t(x^2))
				y=H%*%y2
				r=c(1:d)
				r=sqrt(norm(d,4)/(2*(r+1)*(r+2)))
				y=y*r
				X=t(t(x))%*%t(y)
				for (i in 1:d){
					E[ , ,i]=X*x[i]
				}
				R=E[A]
				return(R)
}

f15=function(x,d,A,H){ 	y=t(x^2)
				y=H%*%t(y)
				y=y*x
				r=c(1:d)
				r=sqrt(norm(d,4)/(6*(r-1)*(r+2)))
				y=y*r
				X=y%*%x
				R=X[A]
				return(R)
}
