reset
set output "Grafica.ps"
a=76.8
b=-17.3
c=-490
f(x)=a+b*x+c*x**2
fit f(x) "Caida.dat" using 1:2 via a,b,c
update "Parametros.txt"

set multiplot
set yrange [0:80]
plot[0:0.35] "Caida.dat" using 1:2
plot[0:0.35] f(x)