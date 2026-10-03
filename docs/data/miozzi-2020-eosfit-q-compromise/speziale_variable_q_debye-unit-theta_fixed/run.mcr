log /Users/clemens/Programming/peritheos/docs/data/miozzi-2020-eosfit-q-compromise/speziale_variable_q_debye-unit-theta_fixed/run.log
read /Users/clemens/Programming/peritheos/docs/data/miozzi-2020-eosfit-q-compromise/speziale_variable_q_debye-unit-theta_fixed/cold.dat
input
pr
2
3
6.8682515367799999
129
6.24
pscale
GPa
vscale
cm^3/mol
x
fit
n
y
y
y
n
n

y
n
save /Users/clemens/Programming/peritheos/docs/data/miozzi-2020-eosfit-q-compromise/speziale_variable_q_debye-unit-theta_fixed/stage1.eos
y
clear
y
read /Users/clemens/Programming/peritheos/docs/data/miozzi-2020-eosfit-q-compromise/speziale_variable_q_debye-unit-theta_fixed/all.dat
input
th
7
300
y

420
1
1.11
x
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
save /Users/clemens/Programming/peritheos/docs/data/miozzi-2020-eosfit-q-compromise/speziale_variable_q_debye-unit-theta_fixed/stage2.eos
y
fit
n
y
y
y
n
y
n
n
n

y
n
save /Users/clemens/Programming/peritheos/docs/data/miozzi-2020-eosfit-q-compromise/speziale_variable_q_debye-unit-theta_fixed/stage3.eos
y
