log /Users/clemens/Programming/peritheos/docs/data/fes-eosfit7c-staged-q-compromise/errors/run.log
read /Users/clemens/Programming/peritheos/docs/data/fes-eosfit7c-staged-q-compromise/errors/cold.dat
input
pr
2
3
14.89877624024
148
4.53
pscale
GPa
vscale
cm^3/mol
x
save /Users/clemens/Programming/peritheos/docs/data/fes-eosfit7c-staged-q-compromise/errors/cold-initial.eos
y
list
fit
n
y
y
y
y
y

y
n
save /Users/clemens/Programming/peritheos/docs/data/fes-eosfit7c-staged-q-compromise/errors/cold-fitted.eos
y
clear
y
read /Users/clemens/Programming/peritheos/docs/data/fes-eosfit7c-staged-q-compromise/errors/all.dat
input
th
7
300
y

417
2
2.42
x
save /Users/clemens/Programming/peritheos/docs/data/fes-eosfit7c-staged-q-compromise/errors/thermal-initial.eos
y
list
fit
n
n
n
n
n
y
y
y
y

y
n
save /Users/clemens/Programming/peritheos/docs/data/fes-eosfit7c-staged-q-compromise/errors/fitted.eos
y
