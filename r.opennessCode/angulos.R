asp2dxdy=function(asp){
      if (asp<=3){dy=-1;dx=-1*(asp-2);}
      if (asp==4 | asp==8){dy=0;dx=(asp-6)/2;}
      if (asp>=5 & asp<=7){dy=1;dx=asp-6;}
      return(c(dx,dy))
}

ang=seq(0,340,20)
ope=runif(length(ang))

cosg=pmax(0,cos((ang+da)*3.14159/180))
plot(ang,cosg,xlim=c(0,360))

mc=ms=c()
h=0
for (da in ang){
    h=h+1;
    cosg=pmax(0,cos((ang+da)*3.14159/180))
    sing=pmax(0,-cos((ang+da)*3.14159/180))  
    wc=cosg/sum(cosg)
    ws=sing/sum(sing)
    mc[h]=crossprod(ope,wc)
    ms[h]=crossprod(ope,ws)
}

plot(ang,mc,type="l")
lines(ang,ms,type="l",col="red")
plot(ang,mc/ms)

ang=seq(0,360,by=.1)
a=45
plot(ang,cos(2*pi*(ang-a)/180))
lines(ang,cos(2*pi*(ang-a-90)/180),col="red")
abline(v=45)
abline(v=135)
abline(v=315)
abline(v=225)

difang=function(a1,a2){
  if(abs(a2-a1)<180) d=abs(a2-a1) else d=360-abs(a2-a1)
  return(d)
}

vdif=vectorize(difang)


ang=seq(0,360,by=.5)
a=130
w1=pmax(0,(Vdif(ang,a)<45)*cos(2*pi*Vdif(ang,a)/180))
w2=pmax(0,(Vdif(ang,a)>45)*cos(2*pi*Vdif(ang,a)/180))
w3=pmax(0,cos(2*pi*(ang-a-90)/180))
plot(ang,w1)
lines(ang,w2,col="red")
lines(ang,w3,col="blue")
abline(v=45)
abline(v=135)
abline(v=315)
abline(v=225)






