log /Users/clemens/Programming/peritheos/docs/data/fes-author-staged/reproduced-errors-tr298/run.log
input
load /Users/clemens/Programming/peritheos/docs/data/fes-author-eosfit/sources/FeS6.eos
th
0
x
read /Users/clemens/Programming/peritheos/docs/data/fes-author-staged/reproduced-errors-tr298-inputs/cold.dat
fit
n
y
y
y
y
y
y

y
n
save /Users/clemens/Programming/peritheos/docs/data/fes-author-staged/reproduced-errors-tr298/cold-fitted.eos
y
clear
y
input
th
7
298
n

417
2
2.42
1
x
read /Users/clemens/Programming/peritheos/docs/data/fes-author-staged/reproduced-errors-tr298-inputs/thermal.dat
save /Users/clemens/Programming/peritheos/docs/data/fes-author-staged/reproduced-errors-tr298/thermal-initial.eos
y
fit
n
n
n
n
n
y
n
y
y
y

y
n
save /Users/clemens/Programming/peritheos/docs/data/fes-author-staged/reproduced-errors-tr298/fitted.eos
y
exit
