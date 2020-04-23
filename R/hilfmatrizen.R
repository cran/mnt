
##Hilfsmatrixen f?r die Berechnung von harmonischen Kugelfunktionen



A2=function(d){
	A=array(F,dim=c(d,d))
	for (i in 1:(d-1)) {
		for (j in (i+1):d) {A[i,j]=T
		}
	}
	return(A)
}

A4=function(d){
	A=array(F,dim=c(d,d,d))
	if(d>2){
	for (i in 1:(d-2)) {
		for (j in (i+1):(d-1)) {
			for (k in (j+1):d) { A[i,j,k]=T
			}
		}
	}}
	return(A)
}
A5=function(d){
	A=array(F,dim=c(d,d))
	if (d>2){
	for (i in 2:(d-1)) {
		for (j in (i+1):d) {A[i,j]=T
		}
	}}
	return(A)
}

A13=function(d){
	A=array(F,dim=c(d,d,d))
	if (d>3){
	for (i in 2:(d-2)) {
		for (j in (i+1):(d-1)) {
			for (k in (j+1):d) { A[i,j,k]=T
			}
		}
	}}
	return(A)
}


A12=function(d){
	A=array(F,dim=c(d,d,d,d))
	if(d>3){
	for (i in 1:(d-3)) {
		for (j in (i+1):(d-2)) {
			for (k in (j+1):(d-1)) {
				for (l in (k+1):d){ A[i,j,k,l]=T
				}
			}
		}
	}}
	return(A)
}

H3=function(d){
	I=matrix(1,d,d)
	S=lower.tri(I)
	S=S*I
	diag(S)=c(0:-(d-1))
	return(S)}



H6=function(d){
	I=matrix(1,d,d)
	S=lower.tri(I)
	S=S*I
	diag(S)=c(-2:-(d+1))
	return(S)}

H7=function(d){
	I=matrix(1,d,d)
	S=lower.tri(I)
	S=S*I
	return(S)}

H9=function(d){
	I=matrix(1,d,d)
	S=lower.tri(I)
	S=S*I*3
	diag(S)=c(-2:-(d+1))
	return(S)}

H10=function(d){
	I=matrix(1,d,d)
	S=lower.tri(I)
	S=S*I
	diag(S)=c(-4:-(d+3))
	return(S)}

H15=function(d){
	I=matrix(1,d,d)
	S=lower.tri(I)
	S=S*I*3
	diag(S)=c(0:-(d-1))
	return(S)}
