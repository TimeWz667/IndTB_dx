U = 0
Asym = 1
Sym = 2
ExCS = 3
TxPub = 4
TxPriOnPub = 5
TxPriOnPri = 6

FLat = 7
SLat = 8
# RLowPub = 9
RHighPub = 10
RStPub = 11
# RLowPri = 12
RHighPri = 13
RStPri = 14
# FLatM72 = 15
# SLatM72 = 16
# UBcg = 17
# FLatBcg = 18
# SLatBcg = 19
# FLatBoth = 20
# SLatBoth = 21

N_States = 15

Infectious = [Asym, Sym, ExCS]
PTB = [Asym, Sym, ExCS, TxPub, TxPriOnPub, TxPriOnPri]
UtTB = [Asym, Sym, ExCS]
LTBI = [FLat, SLat, RHighPub, RStPub, RHighPri, RStPri]
# LTBI += [FLatM72, SLatM72, FLatBcg, SLatBcg, FLatBoth, SLatBoth]

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
