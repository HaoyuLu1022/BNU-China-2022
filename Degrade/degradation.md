# Degradation 

## Overview

In degradation module,  the protein related to the degradation has been formed. Protein CYP2C9 and CYP2C19 Δ9-THC is converted to 11-OH under the co catalysis of  Δ9-THC is converted into Δ 9-THC-COOH, finally formed under the catalysis of UGT Δ 9-THC-COO-glu. 

（degradation）



To model the degradation of Δ9-THC, we will use the Michaelis-Menten rate law, which used a quasi-steady-state assumption to simplify the following enzymatic reaction process:
$$
Substrate (S) + Enzyme (E) ⇋ Enzyme-substrate complex ⇀ Product (P) + Enzyme (E)
$$
This was converted into a quicker and more usable form which described the reaction rate of a substrate becoming a product in the following way:
$$
\frac{d[S]}{dt} = (Vmax[S])/(Km+[S])
$$
V_{max} is the maximum reaction rate based on the amount of enzyme present, K_{m} is the amount of substrate needed in order to achieve half of the maximum reaction rate, and [S] is the concentration. We normalized V_{max}, K_{m} and [S] into the same units through dimensional analysis. V_{max} can be converted to units of nmol·min-1·mg-1. K_{m} and [S] are expressed in units of uM (micromolar, or micromoles per liter volume).



Parameter tables are provided as followed, which state the Vmax and Km values as they were found in the literature. （table）

| Pathway                                     | Enzyme(s)  | Vm               | Km       |           |
| ------------------------------------------- | ---------- | ---------------- | -------- | --------- |
| 11-OH-THC                                   | CYP2C9     | 624 pmol/min/mg  | 0.07 uM  | [1]       |
| 11-OH-THC                                   | CYP2C9     | 6.49[nd]         | 2uM      | [2]       |
| 11-OH-THC                                   | CYP2C9     | 1.1 nmol/min/mg  | 70 nM    | [3]       |
| COOH-THC                                    | CYP2C9     | 54 pmol/min/mg   | 0.50 uM  | [1]       |
| 11-OH-THC metabolite formation(THC-COO-GLU) | UGT2B7/1A9 | 343 pmol/min/mg  | 0.64 uM  | [1]       |
| COOH-THC to THC-COO-glu                     | UGT1A1     | 0.22 nmol/min/mg | 170 uM   | [4]       |
| COOH-THC to THC-COO-glu                     | UGT1A3     | 0.68nmol/min/mg  | 68 uM    | [4]       |
| COOH-THC to THC-COO-glu                     | UGT1A1     | 1.10 nmol/min/mg | 118.3 uM | [5] 60min |
| COOH-THC to THC-COO-glu                     | UGT1A3     | 2.27nmol/min/mg  | 77.1 uM  | [5] 60min |



These steps can be represented in the following equations:


$$
\begin{align}
&\frac{d[\varDelta9-THC]}{dt}=-\frac{V_{m1}[\varDelta9-THC]}{K_{m1}+[\varDelta9-THC]}\\
&\frac{d[11-OH-\varDelta9-THC]}{dt}=\frac{V_{m1}[\varDelta9-THC]}{K_{m1}+[\varDelta9-THC]}-\frac{V_{m2}[11-OH-\varDelta9-THC]}{K_{m1}+[11-OH-\varDelta9-THC]}\\
&\frac{d[\varDelta9-THC-COOH]}{dt}=\frac{V_{m2}[11-OH-\varDelta9-THC]}{K_{m1}+[11-OH-\varDelta9-THC]}-\frac{V_{m3}[\varDelta9-THC-COOH]}{K_{m3}+[\varDelta9-THC-COOH]}\\
&\frac{d[THC-\varDelta9-COO-glu]}{dt}=\frac{V_{m3}[\varDelta9-THC-COOH]}{K_{m3}+[\varDelta9-THC-COOH]}
\end{align}
$$



## Result





## References

[1]:Quantifying Hepatic Enzyme Kinetics of (-)-∆9-Tetrahydrocannabinol (THC) and Its Psychoactive Metabolite, 11-OH-THC, through In Vitro Modeling

[2]:CYP2C-catalyzed delta(9)-tetrahydrocannabinol metabolism: Kinetics, pharmacogenetics and interaction with phenytoin

[3]:Characterizing and Quantifying Extrahepatic Metabolism of (-)-Δ 9-Tetrahydrocannabinol (THC) and Its Psychoactive Metabolite, (±)-11-Hydroxy-Δ 9-THC (11-OH-THC)

[4]Characterization of Human Hepatic and Extrahepatic UDP-Glucuronosyltransferase Enzymes Involved in the Metabolism of Classic Cannabinoids

[5]:Impact of oral Cannabis on driving skills and genetic vulnerability to psychotic symptoms