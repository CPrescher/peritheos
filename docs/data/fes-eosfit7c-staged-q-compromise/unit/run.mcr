log /Users/clemens/Programming/peritheos/docs/data/fes-eosfit7c-staged-q-compromise/unit/run.log
read /Users/clemens/Programming/peritheos/docs/data/fes-eosfit7c-staged-q-compromise/unit/cold.dat
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
save /Users/clemens/Programming/peritheos/docs/data/fes-eosfit7c-staged-q-compromise/unit/cold-initial.eos
y
list
fit
n
y
y
y
n
n

y
n
save /Users/clemens/Programming/peritheos/docs/data/fes-eosfit7c-staged-q-compromise/unit/cold-fitted.eos
y
clear
y
read /Users/clemens/Programming/peritheos/docs/data/fes-eosfit7c-staged-q-compromise/unit/all.dat
input
th
7
300
y

417
2
2.42
x
save /Users/clemens/Programming/peritheos/docs/data/fes-eosfit7c-staged-q-compromise/unit/thermal-initial.eos
y
list
fit
n
n
n
n
n
y
n
n
n

y
n
save /Users/clemens/Programming/peritheos/docs/data/fes-eosfit7c-staged-q-compromise/unit/fitted.eos
y
