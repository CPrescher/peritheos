log /Users/clemens/Programming/peritheos/docs/data/fes-author-staged/published-errors-tr300/run.log
input
load /Users/clemens/Programming/peritheos/docs/data/fes-author-eosfit/sources/FeS6.eos
th
0
pr
2
3
15.4
115.5
4.99
x
read /Users/clemens/Programming/peritheos/docs/data/fes-author-staged/published-errors-tr300-inputs/cold.dat
save /Users/clemens/Programming/peritheos/docs/data/fes-author-staged/published-errors-tr300/cold-fitted.eos
y
clear
y
input
th
7
300
n

417
2
2.42
1
x
read /Users/clemens/Programming/peritheos/docs/data/fes-author-staged/published-errors-tr300-inputs/thermal.dat
save /Users/clemens/Programming/peritheos/docs/data/fes-author-staged/published-errors-tr300/thermal-initial.eos
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
save /Users/clemens/Programming/peritheos/docs/data/fes-author-staged/published-errors-tr300/fitted.eos
y
exit
