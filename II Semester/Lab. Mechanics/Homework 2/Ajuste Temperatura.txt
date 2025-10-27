reset
set output "Grafica.png"
a=-10
b=70
f(x)=a*x+b
fit f(x) "Temperatura.dat" using 1:2 via a,b
update "Parametros(T).txt"

set multiplot
set yrange [0:75]
plot[0:7.5] "Temperatura.dat" using 1:2
plot[0:7.5] f(x)