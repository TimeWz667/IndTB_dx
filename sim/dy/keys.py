U = 0
Asym = 1
Sym = 2
ExCS = 3
ReCS = 4
TxPub = 5
TxPriOnPub = 6
TxPriOnPri = 7

TxNewPub = 8
TxNewPriOnPub = 9
TxNewPriOnPri = 10

FLat = 11
SLat = 12
RHighPub = 13
RStPub = 14
RHighPri = 15
RStPri = 16

N_States = 17

Infectious = [Asym, Sym, ExCS, ReCS, TxNewPub, TxNewPriOnPub, TxNewPriOnPri]
PTB = Infectious + [TxPub, TxPriOnPub, TxPriOnPri]
UtTB = [Asym, Sym, ExCS]
LTBI = [FLat, SLat, RHighPub, RStPub, RHighPri, RStPri]

DS = 0
DR = 1
FR = 2
N_States_R = 3

N_States_A = 8

A_Inc = 0
A_IncRecent = 1
A_IncRemote = 2
A_IncTreated = 3
A_IncTreatedPub = 4
A_Mor = 5
A_NotiPub = 6
A_NotiPri = 7
A_ACF = 8

N_Aux = 9
