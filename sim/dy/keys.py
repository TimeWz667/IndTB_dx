U = 0
Asym = 1
Sym = 2
ExCS = 3
ReCS = 4
TxPub = 5
TxPriOnPub = 6
TxPriOnPri = 7
TxsPub = 8
TxsPriOnPub = 9


FLat = 10
SLat = 11
RHighPub = 12
RStPub = 13
RHighPri = 14
RStPri = 15

N_States = 16

Infectious = [Asym, Sym, ExCS, ReCS]
PTB = Infectious + [TxPub, TxPriOnPub, TxPriOnPri, TxsPub, TxsPriOnPub]
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
