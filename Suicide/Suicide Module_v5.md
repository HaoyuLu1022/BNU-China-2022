# Suicide Module

## Overview

In order to guarantee biosafety of the project, 2 suicide pathways, active and passive respectively, are proposed and implemented in our project.

First, the arabinose-induced suicide pathway is an ODE-based simulation model incorporating spatial aspects, with which we describe the process of arabinose solution diffusing, ccdB expression triggering and cell deaths inside the intestine, and further predict the suicide effect *in vivo*. 

And for the second model, we adopt ODE to simulate and calculate the time all of our bacteria dies, attempting to validate that the suicide decided by leakage to the outside environment happens not long after, thus preventing any possible contamination to it.

## Arabinose-induced Suicide Pathway

### Assumptions

-   We assume that diffusion happens on an ideal surface and therefore neglect the geometry features of the rugged intestinal wall, so that the diffusion process can even out the concentration everywhere.
-   Loss of L-arabinose in the solution when taken in is hard to be exactly quantified since human digestive system is too complicated to simulate. Thus we only consider the decent concentration of L-arabinose solution directly inside intestine, and leave out the reduction or dilution effect in the whole digestive system.

### Description

First, there are chances that our users may want to eliminate the bacteria implanted in their intestine, i.e., to **"actively"** remove it. Hence, we make our bacteria sensitive to the harmless L-arabinose ingestion, which we humans have a low content of inside our body and is scarce in our everyday diet so the killing switch won't be activated accidentally. Specifically, we put the L-arabinose inducible promoter, araC gene, and the gene encoding DNA replication inhibitor (ccdB) into our pathway. Once appropriate amount of L-arabinose is taken in, it will bind with the araC protein produced by our bacteria, and activate the downstream expression of ccdB. After some time, the bacteria will be dead due to the rising amount of ccdB they themselves generate.
$$
\begin{aligned}
&\longrightarrow mRNA_{araC}\\
mRNA_{araC}&\longrightarrow araC\\
A+araC&\longrightarrow A\cdot araC\\
araC &\longrightarrow \empty\\
A\cdot araC&\longrightarrow A\cdot araC +mRNA_{ccdB}\\
mRNA_{ccdB}&\longrightarrow ccdB\\
ccdB&\longrightarrow \empty
\end{aligned}
$$
Considering that L-arabinose solution won't spread evenly in our intestine, we introduce a **reaction-diffusion model** to simulate the process of L-arabinose diffusion, ccdB generation and cell death. We conduct the ODE-based simulation in each grid over a 2-dimensional grid map to represent the intestinal environment based on **Fick's 2nd Law**. Since the distribution of arabinose solution is uneven during its diffusion process, the level of araC at different spots will be different of course, and therefore the time each cell is triggered to express ccdB and dies subsequently is different as well. 

Moreover, according to work from Edwards H et al. [1], araC expression follows **Hill kinetics**, and we refer to part of their kinetics model in our model. To be more specific, L-arabinose taken in binds to the protein araC which used to be a repressor and then activates the downstream gene expression, i.e., pBAD promoter; concentration of araC declines as L-arabinose binds to it; the expression of ccdB is activated by the upstream binding of L-arabinose and araC. We model these 3 biochemistry processes as **hill functions** to offer more precise depiction of our arabinose-induced suicide pathway.

$$
\frac{d\left[A\right]}{dt}=D_A\nabla^2[A]-\frac{k_A\left[araC\right]\left[A\right]^m}{K_A^m+\left[A\right]^m}\\
\frac{d\left[araC\right]}{dt}=\frac{\alpha_{araC}}{\left(\frac{\left[A\right]}{K_{RA}}\right)^p+1}-d_{araC}\left[araC\right]\\
\frac{d[mRNA_{ccdB}]}{dt}=\frac{\gamma_{ccdB}C_N}{\left(\frac{\left[araC\right]}{K_{RaraC}}\right)^q+1}\times flag-d_{mc}[mRNA_{ccdB}]\\
\frac{d\left[ccdB\right]}{dt}=\beta_{ccdB}[mRNA_{ccdB}]-d_{ccdB}\left[ccdB\right]\\
$$
where $flag = 1$ if $\frac{k_A\left[araC\right]\left[A\right]^m}{K_A^m+\left[A\right]^m}>0$ and $=0$ otherwise, indicating that ccdB expression is triggered only when there are araC-arabinose complex forming in the system. The rest of the parameters are stated in the table below.

| Name             | Description                                                 | Value           | Unit          | Reference     |
| ---------------- | ----------------------------------------------------------- | --------------- | ------------- | ------------- |
| $k_A$            | Rate constants for reactant A                               | 0.5             |               | [1]           |
| $K_A$            | Equilibrium constants for reactant A                        | 10              |               | [1]           |
| $\alpha_{araC}$  | Expression rate of protein araC                             | 1               | $s^{-1}$      | [1]           |
| $\gamma _{ccdB}$ | Transcription rate of ccdB                                  | 0.1333          | $s^{-1}$      | [4]           |
| $C_N$            | plasmid number of ccdB                                      | 600             |               | experiment    |
| $d_{mc}$         | Degradation rate of $mRNA_{ccdB}$                           | 0.00203         | $s^{-1}$      | [4]           |
| $\beta_{ccdB}$   | Translation rate of ccdB                                    | 0.033           | $s^{-1}$      | [4]           |
| $K_{RA}$         | Repression equilibrium constant due to repression from A    | 0.2             |               | [1]           |
| $k_{Rarac}$      | Repression equilibrium constant due to repression from araC | 0.1             |               | [1]           |
| $d_{araC}$       | Degradation rate of araC                                    | 0.1             |               | [1]           |
| $d_{ccdB}$       | Degradation rate of ccdB                                    | 0.00115524      | $s^{-1}$      | [4]           |
| $A_0$            | Initial concentration of arabinose                          | 6.66e-6         | $mmol/mm^{3}$ | experiment    |
| $X_0$            | Initial number of *E.coil*                                  | $100\times 100$ |               | estimated     |
| $m$              | Hill coefficient                                            | 2               |               | [1]           |
| $p$              | Hill coefficient                                            | 2               |               | [1]           |
| $q$              | Hill coefficient                                            | 1               |               | [1]           |
| $D$              | diffusion rate for L-arabinose solution                     | 1.037e-3        | $mm^2/s$      | [2]/estimated |

It needs to be noted that the mechanism for ccdB to kill cells is rather complicated for any possible quantification. So we assume that when concentration of ccdB reaches a certain threshold, cells would die in some period of time. And after repeated discussion and lots of attempts to fit the data gained from experiment, we are able to make the reasonable estimation about these 2 important parameters—threshold for ccdB, and the time cells die after ccdB reaching that threshold—that they are 1e-3 mol/mL and 1800s respectively.

### Results

<img src="https://github.com/HalveLuve/Images/blob/master/uPic/arabinose1.png?raw=true" style="zoom: 50%;" >

<img src="https://github.com/HalveLuve/Images/blob/master/uPic/arabinose2.png?raw=true" style="zoom: 67%;" >

In our 100x100 $mm^2$ grid, 5000s simulation, the number of live cells at first increases at a slow rate until around 1800s, and then drops greatly in the later 1800s of the first hour, which is in line with our experiment results and implies that the arabinose-induced cell suicide would be very responsive and quick *in vivo*.

Due to the restrictions of our computing performance, we can't go any further than 5000s in time or 100x100 $mm^2$ in space, otherwise we will run out of memory and the simulation will cost too much time as well. But the current result has already depicted a promising future of our arabinose-induced suicide pathway design.

## mazEF-mediated Suicide Pathway

### Assumptions

-   Unlike the diffusion process of solutions, thermal conduction is very quickly, and hence assume that once leakage occurs, every cell will sense the sudden change in temperature immediately and simultaneously.
-   Regulation at the translation level can be ignored.

### Description

#### Simulation Model

Despite the possibility that our users want to devitalize the intestinal bacteria "actively", it is fact that they may leak into the external environment, causing contamination. So we will have to make the bacteria able to **"passively"** dead when such leakage happens. And this is where our $mazEF$ mediation comes in. Generally, it's a toxin-antitoxin mediation system based on temperature. And in this system, $mazF$ is a toxic protein that can cleave a certain sequence of $mRNA$, and then triggers programmed cell death, while $MazE$ act as an antitoxin that can bind to $mazF$ and inhibit its function. And we suppose that once the concentration of $maxF$ reaches a certain threshold, the cell death effect will be irreversible. In addition, an RBS is introduced into this system to control the expression of $mazE$ based on temperature. The reaction system is shown below.
$$
\begin{aligned}
&\xrightarrow[]{mazF} mRNA_{mazF}\\
&\xrightarrow[]{mazF} mRNA_{mazE}\\
mRNA_{mazF}&\longrightarrow mazF\\
mRNA_{mazE}&\longrightarrow mazE\\
mazE+mazF&\underset{d_{a2}}{\stackrel{a}{\rightleftharpoons}} mazEF\\
mRNA_{mazF}&\xrightarrow[]{d_m} \empty\\
mRNA_{mazE}&\xrightarrow[]{d_m} \empty\\
mazF&\xrightarrow[]{d_c} \empty\\
mazE&\xrightarrow[]{d_a} \empty\\
mazEF&\xrightarrow[]{d_c} \empty\\
\end{aligned}
$$
And we can derive the following equations referring to XJTU-China's work in 2021.
$$
\frac{d\left[mRNA_{MazF}\right]}{dt}=Dk_M\varphi_1\left(T\right)-\frac{d_{large}\left(\beta\left[MazF\right]\right)^2}{\left(\beta\left[MazF\right]\right)^2+K_t^2}\left[mRNA_{MazF}\right]-d_m\left[mRNA_{MazF}\right]\\
\frac{d\left[mRNA_{MazE}\right]}{dt}=Dk_M\varphi_2\left(T\right)-\frac{d_{large}\left(\beta\left[MazF\right]\right)^2}{\left(\beta\left[MazF\right]\right)^2+K_t^2}\left[mRNA_{MazE}\right]-d_m\left[mRNA_{MazE}\right]\\
\frac{d\left[MazF\right]}{dt}=b_1\left[mRNA_{MazF}\right]-a\left[MazF\right]\left[MazE\right]-d_c\left[MazF\right]+d_{a2}\left[MazEF\right]\\
\frac{d\left[MazE\right]}{dt}=b_2\left[mRNA_{MazE}\right]-a\left[MazF\right]\left[MazE\right]-d_a\left[MazE\right]\\
\frac{d\left[MazEF\right]}{dt}=a\left[MazF\right]\left[MazE\right]-d_c\left[MazEF\right]-d_{a2}\left[MazEF\right]
$$
where $\varphi_1(T)$ and $\varphi_2(T)$ are the expression ratio (RBS efficiency) related to temperature for $mazF$ and $mazE$ respectively, and need to be calculated based on experiment result. 

-   For the former, $mazF$ is rather stable in terms of changes in temperature, and there is nothing that will actively inhibit its expression as well, like RBS—but it doesn't mean that its expression keeps a constant level, and one of the reasons is that activity of related enzymes are still subject to temperature fluctuations. 
-   For the latter, however, $mazE$ is more labile, i.e., its degradation rate is much higher than $mazF$ overall. Also, the RBS we add to the system will slow down or even stop its expression when it senses that temperature drops below the normal body temperature. 

This is why $\varphi_1(T)$ and $\varphi_2(T)$ are introduced in our simulation and needed to be fitted. And as for $K_t$, team XJTU-China adopted 15 as the threshold for cleavage of $mRNA$ in their 2020 Model part, which may suit their case well but is pretty likely to work improperly in our situation. So it has to be estimated.

The rest of the parameters are listed below.

| Name        | Description                                   | Value           | Unit     | Reference  |
| ----------- | --------------------------------------------- | --------------- | -------- | ---------- |
| $D$         | Concentration of the promoter                 | 3049            |          | experiment |
| $k_M$       | Maximum initiation rate of the promoter       | 0.03            | $s^{-1}$ | [3]        |
| $b_1$       | $mazF$ translation rate                       | 0.009           | $s^{-1}$ | [3]        |
| $b_2$       | mazE translation rate                         | 0.0317          | $s^{-1}$ | estimated  |
| $a$         | Binding rate of mazE and $mazF$               | 0.000365        | $s^{-1}$ | [3]        |
| $d_m$       | Degradation rate of mRNA                      | 0.002           | $s^{-1}$ | [3]        |
| $d_c$       | Degradation rate due to cell division         | 0.00028         | $s^{-1}$ | [3]        |
| $d_a$       | Degradation rate of mazE                      | 0.00231         | $s^{-1}$ | [3]        |
| $d_{a2}$    | Degradation rate of mazE within $mazEF$       | 0.000231        | $s^{-1}$ | [3]        |
| $d_{large}$ | Increased mRNA degradation rate due to $mazF$ | 0.198           | $s^{-1}$ | [3]        |
| $\beta$     | Concentration factor of $mazF$                | 0.1             | $s^{-1}$ | [3]        |
| $K_t$       | Threshold for cleavage of mRNA                | discussed later |          | estimated  |

However, $mazF$ kills cells, which brings more uncertainty to the current system, thus making it difficult to obtain these values through fitting. So we construct another new plasmid, using $GFP$ to replace the $mazF$ segment, and hence we can fit and obtain $\varphi_1(T)$ and $\varphi_2(T)$ by establishing a modified model where there is no need to consider cell deaths.

Basically, what we need to do is calculate $GFP$ concentration/amount based on the strength of fluorescence detected. Therefore, the modified model contains 2 parts: first, $GFP$ expression needs to be quantified, and second, the mapping/relation between fluorescence intensity and concentration of $GFP$ needs to be specified.

#### Sub-model 1: GFP Expression

To simulate the $GFP$ expression level, we refer to Stögbauer T et al.'s work [5] and derive the reactions and equations below. 
$$
\begin{aligned}
&\xrightarrow[]{TsR}mRNA_{GFP}\\
mRNA_{GFP}&\xrightarrow[]{TlR}GFP\\
mRNA_{GFP}&\xrightarrow[]{d_{mRNA}}\empty\\
GFP&\xrightarrow[]{k_{mat}}\empty\\
\end{aligned}
$$

$$
\frac{d[mRNA]}{dt}=\frac{k_{ts}\times [TsR]\times [DNA]}{K_s+[DNA]}-d_{mRNA}\times [mRNA]\\
\frac{d[GFP]}{dt}=\frac{k_{tl}\times [TlR]\times [mRNA]}{K_l+[mRNA]}-d_{GFP}\times [GFP]\\
$$

where $k_{ts}$ and $k_{tl}$ represent the rates of transcription and translation respectively, and $d_{mRNA}$ refers to the degradation rate of mRNA for $GFP$. $d_{GFP}$ is the rate at which $GFP$ degrades. It needs to be pointed out that, we define $TsR$ and $TlR$ as all the transcription and translation resources including polymerases, ribosomes, tRNAs, etc., according to Stögbauer T et al.'s work [5] and that we follow their idea to set their initial concentrations to $1nM(nmol/L)$. $K_s$ and $K_l$ are Michaelis constants for $TsR$ and $TlR$ respectively. In addition, Stögbauer T et al. conducted their research *in vitro*, where the transcription and translation resources are limited and consumed but not produced or updated. So in our sub-model we assume that there are always stable contents of $TsR$ and $TlR$ in each single cell.
| Name       | Description                | Value                                   | Unit       | reference  |
| ---------- | -------------------------- | --------------------------------------- | ---------- | ---------- |
| $k_{ts}$   | Transcription rate         | $18.2$                                  | $nM/min$   | [5]        |
| $k_{tl}$   | Translation rate           | $16.1$                                  | $nM/min$   | [5]        |
| $d_{GFP}$  | $GFP$ degradation rate     | $0.156$                                 | $min^{-1}$ | [6]        |
| $d_{mRNA}$ | mRNA degradation rate      | $8.4e-2$                                | $min^{-1}$ | [7]        |
| $K_s$      | Michaelis constant         | $0.8\sim4.2$                            | $nM$       | [8]        |
| $K_l$      | Michaelis constant         | $14\sim30$（previous work, *in vitro*） | $nM$       | [9]        |
| $[DNA]$    | Template DNA concentration | $998.14\sim 9981.4$                     | $nM$       | experiment |
| $[TsR]$    | Transcription resource     | $1$                                     |            | [5]        |
| $[TlR]$    | Translation resource       | $1$                                     |            | [5]        |

#### Sub-model 2: Fluorescence intensity and concentration of GFP

Since the data collected is fluorescence intensity, and concentration of $GFP$ is difficult to obtain by experiment, we need to establish another sub-model to describe the mapping/relation between them. We refer to the principle of Quantitative Fluorescence Analysis and gain a much deeper understanding about our problem. 

Basically, the method relies on the formula
$$
\begin{aligned}
F&=KI_0(1-10^{-ECl})\\
&=KI_0(1-e^{-2.3\times ECl})
\end{aligned}
$$
where $K$ is a constant that depends on fluorescence efficiency, $I_0$ means the initial intensity of the fluorescence, and $ECl$ together calculates the absorbance based on Lambert-Beer's Law, among which $E$ refers to the extinction coefficient, $C$ is the concentration of $GFP$ and $l$ stands for the optical length. Then we transform the equations into forms shown below through power series expansion.
$$
\begin{aligned}
F&=KI_0(1-\sum\limits_{i=0}\frac{(-2.3\times ECl)^i}{i!})\\
&=KI_0(2.3\times ECl+\sum\limits_{i=2}\frac{(-2.3\times ECl)^i}{i!})
\end{aligned}
$$
Now consider the term $ECl$: $E$ and $l$ are constants that are irrelevant to concentration, while $C$, $GFP$'s concentration, is the only variable in this equation. When $C$ is low enough to make $ECl < 0.05$, then it's safe to eliminate the term $\sum\limits_{i=2}\frac{(-2.3\times ECl)^i}{i!}$ from the equation while we won't make it inaccurate, thus only keeping the 1-order term $2.3\times ECl$. 

So far, we have derived the linear approximation form of this formula, and for further simplification, we replace the rather complex term $2.3KI_0El$ with $K^{'}$ since they are all constants and now we have the equation
$$
F\approx 2.3KI_0ECl=K^{'}C
$$
Therefore, we have proven that when concentration of $GFP$ is low enough ($ECl < 0.05$ in this case), we can draw the safe conclusion that the fluorescence intensity is linear with the concentration of $GFP$. Furthermore, this conclusion implies that we can simply introduce the **normalized fluorescence intensity** to represent the $GFP$ concentration in our model.

#### Fitting Results

<img src="https://raw.githubusercontent.com/HalveLuve/Images/master/PicGo/image-20221010092929067.png" alt="image-20221010092929067" style="zoom:50%;" />

Through some curve fitting techniques, we managed to find the optimal parameters for $\varphi_1(T)$ and $\varphi_2(T)$, which are 1 and 0.132 respectively. And based on this fitting result, we can further estimate the threshold for cleavage of $mRNA$, namely $K_t$, to be at least above the concentration of $mazF$ when temperature is greater than 37 degrees Celsius so that when inside the human intestinal environment, the bacteria won't be falsely triggered to suicide. And after repeated testing and validation, we eventually estimate it to be at least 90.

### Results

![](https://raw.githubusercontent.com/HalveLuve/Images/master/PicGo/image-20221009010313587.png)

On the left side is the situation where temperature >= 37 ℃, and concentration of $mazF$, the toxin, remains at single digit throughout the 2-hour simulation, which is way much lower than the cleavage threshold. While the simulation result of <= 37 ℃ on the right side implies that $mazF$ can be produced at a pretty fast rate since its degradation rate is slow enough to be ignored, and there is not much $mazE$ that binds to it and inhibits its toxin effect. In this case, the concentration of the antitoxin $mazE$ peaks around the first few tens of seconds and then declines slowly until 2 hours later. Hence it's safe to draw the conclusion that the bacteria can be eliminated immediately when leakage happens.

## References

[1]: Edwards H, Xu P. Unstructured kinetic models to simulate an arabinose switch that decouples cell growth from metabolite production. Synth Syst Biotechnol. 2020 Jul 14;5(3):222-229. doi: 10.1016/j.synbio.2020.07.003. PMID: 32695893; PMCID: PMC7364165.

[2]: Mogi, Noriko et al. “Infinite dilution binary diffusion coefficients for six sugars at 0.1 MPa and temperatures from (273.2 to 353.2) K.” *Journal of Chemical & Engineering Data* 52 (2007): 40-43.Infinite Dilution Binary Diffusion Coefficients for Six Sugars at 0.1 MPa and Temperatures from (273.2 to 353.2) K

[3]: 2020.igem.org/ XJTU-China/Suicide-Model

[4]: Gelens L, Hill L, Vandervelde A, Danckaert J, Loris R (2013) A General Model for Toxin-Antitoxin Module Dynamics Can Explain Persister Cell Formation in *E. coli*. PLoS Comput Biol 9(8): e1003190. https://doi.org/10.1371/journal.pcbi.1003190

[5]: Stögbauer T, Windhager L, Zimmer R, Rädler JO. Experiment and mathematical modeling of gene expression dynamics in a cell-free system. Integr Biol (Camb). 2012 May;4(5):494-501.

[6]: J. A. Megerle, G. Fritz, U. Gerland, K. Jung and J. O. Radler, Biophys. J., 2008, 95, 2103–2115.

[7]: E. Karzbrun, J. Shin, R. Bar-Ziv and V. Noireaux, Phys. Rev. Lett., 2011, 106, 1–4.

[8]: M. Maslak and C. T. Martin, Biochemistry, 1993, 32, 4281–4285.

[9]: S. Takahashi, R. Akita, H. Matsuno, H. Furusawa, Y. Shimizu, T. Ueda and Y. Okahata, ChemBioChem, 2008, 9, 870–873.

